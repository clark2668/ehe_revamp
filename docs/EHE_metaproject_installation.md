### EHE metaproject installation instruction

0. First get your python version from cvmfs.
    ` eval `/cvmfs/icecube.opensciencegrid.org/py3-v4.0.1/setup.sh` `

1. Create a folder for your installation, lets call it $I3_BASE. You can set it in your bash by running `export I3_BASE=/your/fav/path`.
    `cd $I3_BASE`

2. Check out combo release and grab a cup of coffee (takes some time)
	`svn co http://code.icecube.wisc.edu/svn/meta-projects/combo/releases/V00-00-03 src`

3. Go to source directory and check out the relevant EHE projects (more might be added in the future)
    `cd src`
    `svn co http://code.icecube.wisc.edu/svn/projects/c2j-icetray/trunk c2j-icetray`
    `svn co http://code.icecube.wisc.edu/svn/projects/juliet/trunk juliet`
    `svn co http://code.icecube.wisc.edu/svn/projects/juliet-interface/trunk juliet-interface`
    `svn co http://code.icecube.wisc.edu/svn/projects/weighting-module/trunk weighting-module`

4. Clean up the juliet directory and grab the newest source files from github
    `cd juliet`
    `rm -r java_lib`
    `git clone https://github.com/ShigeruYoshida/JULIeT.git java_lib`

5. Go back to $I3_BASE and make a build directory
    `cd $I3_BASE`
    `mkdir build && cd build`

6. Set your $JAVA_HOME and run cmake
	`export JAVA_HOME=/cvmfs/icecube.opensciencegrid.org/users/mmeier/jdk1.7.0_80`
	`cmake -DBUILD_JULIET:BOOL=ON -DI3_EXTRA_TOOLS=jni ../src`


7. Run make. You can use -jN to compile with multiple (N) cores! (Aaand you can grab another coffee!)
    `make`

8. Make sure your juliet can find the proagation matrices. You can adjust the path to the propagation matrix if you're not installing on cobalt. E.g. you can download them from this url: http://www.ppl.phys.chiba-u.jp/JULIeT/
    `cd $I3_BASE/src/juliet/java_lib`
    `mkdir data`
    `ln -s /data/user/mmeier/propMtx/ data/`
    (`ln -s /home/mmeier/data/propMtx data/` should work on grappa)

9. To test juliet and the weighting-module enter an env-shell.sh and run the script test_juliet.py
    `cd $I3_BASE/build`
    `./env-shell.sh`
    On cobalt:
    `python /data/user/mmeier/table_based_sim/juliet_scripts/test_juliet.py --flavor mu --filename /path/to/the/output/file.i3`
    On grappa:
    `python /home/mmeier/data/software/test_ehe/test_juliet.py --flavor mu --filename /path/to/the/output/file.i3`

    Flavor can be (mu, tau, nue, numu, nutau)

