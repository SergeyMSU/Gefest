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


include "elliptic_integral.f90"
include "Solvers.f90"
include "cgod_3D.f90"
include "Storage.f90"
include "MK.f90"
include "Tab1.f90"
include "My_func.f90"
include "Distfunc.f90"



program Gefest
	USE STORAGE
	USE Distfunc
	USE MK
	USE Tab1
	USE Emath
	implicit none
	real(8) :: TT, r1, All_TT
	integer :: i, now
	
	now = 1
	! Variables

	! Body of Gefest
	print *, 'Init'
	call Init_Func(f1)
	call Init_Func(f2)
	call Init_Func(f3)
	call Init_FuncGD(gd1)
	call Tab1_Set(f1)
	call MK_Read_Sig()
	print *, 'Start'
	!call Test_int()

	!print*, "time = ", gd1%time_step
	gd1%time_step = 0.0000001_8
	TT = 0.0
	do while(TT < 0.0) !0.0612
		r1 = gd1%time_step
		gd1%time_step = 100000000.0
		call Start_GD(r1, gd1, 2)
		TT = TT + r1
		print*, "time = ", r1, TT
		r1 = gd1%time_step
		gd1%time_step = 100000000.0
		call Start_GD(r1, gd1, 1)
		TT = TT + r1
		print*, "time = ", r1, TT
	end do

	call Print_GD(gd1, "00001", TT)

	!call Calc_Q(f1, gd1, 1)
	TT = 0.1_8
	All_TT = 0.0_8
	! call Print_fx(f1, -0.5_8, "00003")
	! call Print_fx(f1, 0.5_8, "00001")
	! call Print_rho(f1, "00001")
	! call Start(TT)


	do i = 1, 2
		print*, "________________________ global step = ", i 
		TT = 0.02_8
		call Integrate_Protiv_potoka(f1, f2, f3, TT, now)
		All_TT = All_TT + TT
		print*, "Time = ", TT, All_TT
		call Print_fx(f1, 0.01_8, "00001", All_TT, i)
		call Print_fx(f1, 0.5_8, "00002", All_TT, i)
		call Print_fx(f1, 1.0_8, "00003", All_TT, i)
		call Print_fx(f1, 1.5_8, "00004", All_TT, i)
		call Print_fx(f1, 2.0_8, "00005", All_TT, i)
		call Print_fx(f1, 3.0_8, "00006", All_TT, i)
		call Print_fx(f1, -0.1_8, "00007", All_TT, i)
		call Print_fx(f1, -0.2_8, "00008", All_TT, i)
		call Print_fx(f1, -0.3_8, "00009", All_TT, i)
		call Print_fx(f1, 0.2_8, "00010", All_TT, i)
		call Print_fx(f1, 0.3_8, "00011", All_TT, i)
		call Print_fx(f1, 0.4_8, "00012", All_TT, i)
		call Print_fx(f1, -0.4_8, "00013", All_TT, i)
		call Print_fx(f1, -0.5_8, "00014", All_TT, i)
		call Print_fx(f1, -1.0_8, "00015", All_TT, i)
		call Print_fx(f1, -2.0_8, "00016", All_TT, i)
		call Print_fx(f1, -1.5_8, "00017", All_TT, i)
		call Print_fx(f1, 0.0_8, "00018", All_TT, i)
		! call Print_fx(f1, 1.0_8, "00010")
		! call Print_fx(f1, 2.0_8, "00011")
		! call Print_fx(f1, 3.0_8, "00012")
		! call Print_fx(f1, 4.0_8, "00013")

		! call Print_rho(f1, "00001")
		! call Print_GD(gd1, "00001", All_TT)

		call Save_setka_bin(i)
	end do
	


	! call Print_fx(f1, 0.1_8, "00002")
	! call Print_fx_k(f1, 330, "00002", All_TT)
	! call Print_fx(f1, -0.1_8, "00004")
	! call Print_rho(f1, "00002")


	print*, "par_ch1 = ", par_ch1
	!print*, Tab1_Get(1.0_8, 0.9_8), Tab1_Get(1.0_8, 0.95_8), Tab1_Get(1.0_8, 0.98_8), Tab1_Get(1.0_8, 0.99_8), Tab1_Get(1.0_8, 0.999_8)
	!pause

end program Gefest

