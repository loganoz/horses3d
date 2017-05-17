1.- Execute sh configure.sh
2.- Compile Tests cases with "make all" in each folder 
3.- The binary of each test case is stored in bin folder

----------------------------------------------------------------------------
Changes to master during week 06.03.2017 - 12.03.2017:
	grubio:
		1. STOP 99 in test cases
	arueda:
		1. Time computed as t+dt and not as t_0+k*dt (allows non constant dt)
		2. RK3 residual now computed after time step incrementation
		3. New constructor for TimeIntegrator_t class
		4. Implicit Jacobian-Free Newton-Krylov solver added with new test case under Tests/Euler/JFNK
		5. Expected solutions of test cases changed due to arueda(1) and arueda(2)

----------------------------------------------------------------------------
This modification is added to force an automatic bulding of buildbot and 
test some modifications that have been introduced
