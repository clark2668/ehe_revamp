#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.0/icetray-start
#METAPROJECT: /mnt/home/baclark/IceCube/new_ehe/new/metaproject/build/
from I3Tray import *
from math import *
from os.path import expandvars

from icecube.sim_services.sim_utils.gcd_utils import get_time
from icecube import icetray, dataclasses, dataio
import numpy as np

import os
from optparse import OptionParser

gcd_default = os.path.join(
    '/cvmfs/icecube.opensciencegrid.org/data/GCD',
    'GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz')
filename_default = os.path.join(
    '/data/user/mmeier/table_based_sim/juliet',
    'mu_test_new.i3')

# Setup the parser on top-level so we can get the args when importing
parser = OptionParser()
parser.add_option('-f', '--filename', type='string', dest='filename',
	      default=filename_default)
parser.add_option('-g', '--gcdfile', type='string', dest='gcd_file',
	      default=gcd_default)
# parser.add_option('-s', '--seed', type='int', dest='seed')
parser.add_option('--flavor', dest='flavor',
	      choices=('nue', 'numu', 'nutau', 'mu', 'tau'))
parser.add_option('-r', '--run_number', type='int', default=1,
	      dest='run_number')
parser.add_option('-n', '--n_events', type='int',
	       default=10, dest='n_events')
parser.add_option('--gamma', type='float', default=1.,
	      dest='gamma')
parser.add_option('--energy_min', type='float', default=1e5,
	      dest='e_min')
parser.add_option('--energy_max', type='float', default=1e11,
	      dest='e_max')


load("libc2j-icetray")
load("libicetray")
load("libdataclasses")
load('libdataio')
load("libphys-services")
load("libsim-services")
load("libjuliet-interface")
load("libweighting-module")

rnd_num = 10
JAVA_CLASS_PATH = os.path.expandvars('$I3_SRC/juliet/java_lib/classes')
FLAVOR_TO_PID = {
    'e': 11,
    'nue': 12,
    'mu': 13,
    'numu': 14,
    'tau': 15,
    'nutau': 16
}

FLAVOR_TO_JULIET = {
    'e': 0,
    'nue': 0,
    'mu': 1,
    'numu': 1,
    'tau': 2,
    'nutau': 2
}

MATERIAL = 'ice'  # can also be 'rock'
START_LOCATION = 'cylinder' # can also be 'earth_surface'


def DrivingTime(frame, time):
    if "DrivingTime" in frame:
        del frame["DrivingTime"]
    frame.Put("DrivingTime", time)


def ensure_dir(filename):
    dirname = os.path.dirname(filename)
    if not os.path.isdir(dirname):
        try:
            os.makedirs(dirname, exist_ok=True)
        except:
            print('Directory creation did not work properly, exiting!')
            exit(1)


def main(filename, gcd_file, flavor, run_number,
         n_events, gamma, e_min, e_max):
    # First, make sure the intended output folder exists
    ensure_dir(filename)

    time = get_time(dataio.I3File(gcd_file))
    tray = I3Tray()

    #icetray.logging.set_level(icetray.logging.info)
    icetray.set_log_level(icetray.I3LogLevel.LOG_TRACE)
    #
    # parameters
    #
    # tray.AddService("I3SPRNGRandomServiceFactory","random")(
    #         ("Seed",0),
    #         ("NStreams",500),
    #         ("StreamNum",rnd_num))
    seed = int(FLAVOR_TO_PID[flavor] * 1e5 + run_number)

    tray.AddService('I3GSLRandomServiceFactory', 'random',
		    Seed=seed)

    tray.AddService("I3JavaVMFactory", "java_vm",
		    Options=[expandvars(f"-Djava.class.path={JAVA_CLASS_PATH}"),
			     "-Xms512m", "-Xmx1024m"])

    tray.AddModule("I3InfiniteSource", "somanyevents",
		   Prefix=gcd_file,
		   Stream=icetray.I3Frame.DAQ)

    tray.AddModule(DrivingTime, "dt",
                   time=time,
		   Streams=[icetray.I3Frame.DAQ])

    tray.AddModule("I3MCEventHeaderGenerator", "time-gen",
		   Year=time.utc_year,
		   DAQTime=time.utc_daq_time)

    # Define if the primary is a neutrino (0)
    # or a charged lepton (1) for juliet
    if 'nu' in flavor:
        primary_doublet = 0
        is_neutrino = True
    else:
        primary_doublet = 1
        is_neutrino = False

    # Primary flavor corresponds to the actual lepton flavor
    # e: 0, mu: 1, tau: 2
    primary_flavor = FLAVOR_TO_JULIET[flavor]
    calc_tau_losses = primary_flavor >= 1
    calc_mu_losses = primary_flavor >= 1

    if MATERIAL is 'ice':
        material_id = 0
    elif MATERIAL is 'rock':
        material_id = 1
    else:
        raise ValueError('MATERIAL has to be either ice or rock!')

    if START_LOCATION is 'cylinder':
         start_location_id = 2
    elif START_LOCATION is 'earth_surface':
         start_location_id = 1
    else:
        raise ValueError(
            'START_LOCATION has to be either cylinder or earth_surface!')

    tray.AddModule("I3JulietPrimaryParticleSource", "fillprimary",
		   PrimaryFlavor=primary_flavor,
		   PrimaryDoublet=primary_doublet,
		   EnergyMin=e_min*I3Units.GeV,
		   EnergyMax=e_max*I3Units.GeV,
		   DifferentialPowerLawIndex=gamma,
		   InjectRadius=880.*I3Units.m,
		   NadirMin=0.*I3Units.deg,
		   NadirMax=180.*I3Units.deg,
		   AzimuthMin=0.*I3Units.deg,
		   AzimuthMax=360.*I3Units.deg,
		   NeutrinoInteractionWeight=1)

    tray.AddModule("I3Juliet","propagate",
		   MaterialID=material_id,
		   Flavor=primary_flavor,
		   DoChargedCurrent=is_neutrino,
		   DoNeutralCurrent=is_neutrino,
		   DoMuBremss=calc_mu_losses,
		   DoTauBremss=calc_tau_losses,
		   DoMuKnockOn=0,
		   DoTauKnockOn=0,
		   DoMu2ePairCreation=calc_mu_losses,
		   DoTau2ePairCreation=calc_tau_losses,
		   DoMu2muPairCreation=0,
		   DoTau2muPairCreation=0,
		   DoMu2tauPairCreation=0,
		   DoTau2tauPairCreation=calc_tau_losses,
		   DoMuPhotoNuclear=calc_mu_losses,
		   DoTauPhotoNuclear=calc_tau_losses,
		   DoMuDecay=calc_mu_losses,
		   DoTauDecay=calc_tau_losses,
		   DoGlashowResonance=0,
		   StartLocationID=start_location_id)

    tray.AddModule("I3NullSplitter", "null_split",
		   SubEventStreamName='null_split')

    tray.AddModule('I3WeigherModuleJuliet', 'nu_weight',
		   frameMCWeightName='JulietWeightDict',
		   PrimaryParticleType=0, # 0: neutrino  1: atmosphericMuon           
		   NeutrinoModelList=np.arange(1, 38).tolist(),
		   NDOMthres=0)

    tray.AddModule('I3WeigherModuleJuliet', 'atmo_weight',
		   frameMCWeightName='JulietWeightDictAtmo',
		   frameCosmicRayEnergyWeightName='CosmicRayEnergyDist',
		   PrimaryParticleType=0, # 0: neutrino  1: atmosphericMuon
		   NDOMthres=0,
		   Alpha=1.97, # The Elbert formula's parameter
		   MuonEthreshold=1505, # The Elbert formula's parameter [GeV]
		   CELFlag=False, # use the propagation matrix
		   gzkCutOffFlag=False, # With NO GZK cutoff. But you will fill the cutoff case later anyway
		   AtmFluxWeightName='ElbertModelIC22',
		   widthOfLogEnergyBin=0.05 # width of the logE bin of the CR energy distribution
		   )

    tray.AddModule("I3Writer","writer",
		   streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics,
                            icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
		   filename=filename)

    tray.AddModule("Dump", "dump")
    tray.AddModule("TrashCan", "the can")

    #tray.Execute(103)
    tray.Execute(n_events + 3)
    tray.Finish()


if __name__ == '__main__':
    (options, args) = parser.parse_args()
    filename = options.filename
    gcd_file = options.gcd_file
    # seed = options.seed
    flavor = options.flavor
    run_number = options.run_number
    n_events = options.n_events
    gamma = options.gamma
    e_min = options.e_min
    e_max = options.e_max

    main(filename, gcd_file, flavor, run_number, 
         n_events, gamma, e_min, e_max)

