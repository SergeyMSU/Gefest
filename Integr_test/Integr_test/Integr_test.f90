!  Integr_test.f90 
!
!  FUNCTIONS:
!  Integr_test - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Integr_test
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************
	
MODULE MY
	
	integer :: s1 = 1
	integer :: s2 = 2
	integer :: s3 = 3
	
	contains
	
	subroutine M_K_rand(b)
		real(8), intent(out) :: b
		integer(4):: ic15 = 32768, ic10 = 1024
		integer(4):: mz = 710, my = 17784, mx = 11973
		real(8):: xi = 9.0949470177292824E-13, c = 1.073741824E9
		integer(4) :: i13, i12, i11, ii
		i13 = mz * s1 + my * s2 + mx * s3
		i12 = my * s1 + mx * s2
		i11 = mx * s1
		ii = i11 / ic15
		i12 = i12 + ii
		s1 = i11 - ic15 * ii
		ii = i12 / ic15
		i13 = i13 + ii
		s2 = i12 - ic15 * ii
		s3 = mod(i13,ic10)
		b = xi * (c * s3 + ic15 * s2 + s1)
	end subroutine M_K_rand
	
	function fff(y)
		real(8), intent(in) :: y
		real(8) :: fff
		if(y <= 1.0) fff = 1.0
		if(y > 1.0) fff = 0.0
	end function fff


end MODULE MY

    program Integr_test
    USE MY
    implicit none
	real(8) :: ksi1, ksi2, i, x, y, S

	S = 0.0
    do i = 1, 10000
		call M_K_rand(ksi1)
		x = ksi1
		S = S + x * sin(x)
	end do
	
	S = S/10000

    ! Body of Integr_test
    print *, "1 = ", S, dabs(S - 0.301169)
	
	S = 0.0
    do i = 1, 10000
		call M_K_rand(ksi1)
		x = acos(1.0 - ksi1 * (1.0 - cos(1.0)))
		!print*, x
		!pause
		S = S + x
	end do
	
	S = (1.0 - cos(1.0)) * S/10000
	print *, "2 = ", S, dabs(S - 0.301169)
	
	pause

    end program Integr_test

