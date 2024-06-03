!  Gefest.f90 
!
!  FUNCTIONS:
!  Gefest - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Gefest
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************


include "cgod_3D.f90"
include "Storage.f90"
include "Tab1.f90"
include "My_func.f90"
include "Distfunc.f90"


program Gefest
	USE STORAGE
	USE Distfunc
	USE Tab1
	implicit none
	real(8) :: TT, r1
	integer :: i


	! Variables

	! Body of Gefest
	print *, 'Init', mod(1, 3)
	call Init_Func(f1)
	call Init_Func(f2)
	call Init_Func(f3)
	call Init_FuncGD(gd1)
	call Tab1_Set()
	print *, 'Start'

	!print*, "time = ", gd1%time_step
	gd1%time_step = 0.0000001_8
	TT = 0.0
	do while(TT < 0.0)
		r1 = gd1%time_step
		gd1%time_step = 100000000.0
		call Start_GD(r1, gd1, 2)
		TT = TT + r1
		!print*, "time = ", gd1%time_step
		r1 = gd1%time_step
		gd1%time_step = 100000000.0
		call Start_GD(r1, gd1, 1)
		TT = TT + r1
		!print*, "time = ", gd1%time_step
	end do

	!call Print_GD(gd1, "00001", TT)

	!call Calc_Q(f1, gd1, 1)
	TT = 0.1_8
	! call Print_fx(f1, -0.5_8, "00003")
	! call Print_fx(f1, 0.5_8, "00001")
	! call Print_rho(f1, "00001")
	! call Start(TT)

	call Integrate_Protiv_potoka(f1, f2, f3, TT)



	call Print_fx(f1, 0.1_8, "00002")
	call Print_fx_k(f1, 330, "00002", TT)
	call Print_fx(f1, -0.1_8, "00004")
	call Print_rho(f1, "00002")


	print*, "par_ch1 = ", par_ch1
	!print*, Tab1_Get(1.0_8, 0.9_8), Tab1_Get(1.0_8, 0.95_8), Tab1_Get(1.0_8, 0.98_8), Tab1_Get(1.0_8, 0.99_8), Tab1_Get(1.0_8, 0.999_8)
	pause

end program Gefest

