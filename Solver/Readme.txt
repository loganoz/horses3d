################################################################################
#                                                                              #
#                                                                              #
#            HORSES3D High-Order (DG) Spectral Element Solver                  #
#                                                                              #
#                                                                              #
################################################################################



----------------------------------------------------------------------------
1.- Create a build directory
		$ mkdir build
2.- Set up cmake from build directory
		$ cd build
		$ cmake .. <<options>>
		
		(you can also do this graphically using ccmake ore cmake-gui>>
3.- Install
		$ make install
4.- The program is ready to be used. If you want to test the code, configure
		the test cases and follow the instructions in each test directory.
		$ ./configure
3.- Be sure to add Source/lib to LD_LIBRARY_PATH, so HORSES 3D can be used

----------------------------------------------------------------------------
