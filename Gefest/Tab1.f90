module Tab1
    USE STORAGE
    USE Emath
    USE OMP_LIB
	implicit none 
    !? Затабулированная функция для вычисления интегралла (a - b * cos(phi))

    real(8) :: Tab1_R1 = 20.0_8
    real(8) :: Tab1_R2 = 20.0_8
    integer(4) :: Tab1_n1 = 1000
    integer(4) :: Tab1_n2 = 1000

    real(8), allocatable :: Omega1(:, :)

    contains

    function Tab1_Get(a, b)
        real(8), intent(in) :: a, b
        real(8) :: al1, al2, u1, u2, v1, v2
        real(8) :: Tab1_Get
        integer(4) :: i1, i2, j1, j2
        
        
		i1 = INT(Tab1_n1 * a/Tab1_R1 + 0.5)
		i2 = i1 + 1

        if(i1 < 1) then
            i1 = 1
            i2 = 2
        end if

        if(i2 > Tab1_n1) then
            i2 = Tab1_n1
            i1 = i2 - 1
        end if

        if(a > Tab1_R1) then
            i2 = Tab1_n1
            i1 = i2 - 1
        end if
		
        u1 = Tab1_R1 * (i1 - 0.5)/Tab1_n1
        u2 = Tab1_R1 * (i2 - 0.5)/Tab1_n1


        j1 = INT(Tab1_n2 * b/Tab1_R2 + 0.5)
		j2 = j1 + 1

        if(j1 < 1) then
            j1 = 1
            j2 = 2
        end if

        if(j2 > Tab1_n2) then
            j2 = Tab1_n2
            j1 = j2 - 1
        end if

        if(b > Tab1_R2) then
            j2 = Tab1_n2
            j1 = j2 - 1
        end if
		
        !v1 = Tab1_R2 * (j1 - 0.5)/Tab1_n2
        !v2 = Tab1_R2 * (j2 - 0.5)/Tab1_n2

        !al1 = (a - u1)/(u2 - u1)
        !al2 = (b - v1)/(v2 - v1)


        Tab1_Get = Omega1(i2, j1)

        return
        !if(u1 > v1 .and. u1 > v2) then
		Tab1_Get = (1.0 - al2) * ((1.0 - al1) * Omega1(i1, j1) + al1 * Omega1(i2, j1)) + al2 * ((1.0 - al1) * Omega1(i1, j2) + al1 * Omega1(i2, j2))
        ! else
        !     Tab1_Get = Omega1(i2, j1)
        ! end if
   
    end function Tab1_Get

    subroutine Tab1_Set()

        integer(4) :: i, j, k, n3
        real(8) :: phi, dphi, S, a, b, u, r

        ALLOCATE(Omega1(Tab1_n1, Tab1_n2))
		Omega1 = 0.0

        n3 = 200
        dphi = par_pi/n3

        !$omp parallel

		!$omp do private(j, k, phi, S, a, b, u, r) schedule(dynamic, 3)
        do i = 1, Tab1_n1
            do j = 1, Tab1_n2
                S = 0.0
                a = Tab1_R1 * (i - 0.5)/Tab1_n1
                b = Tab1_R2 * (j - 0.5)/Tab1_n2
                if(a < b) a = b * 1.0000001_8
                
                do k = 1, n3
                    phi = (k - 0.5) * dphi
                    u = sqrt(a - b * cos(phi))
                    S = S + u * sig(u) * dphi
                end do
                Omega1(i, j) = S
                !!r = sqrt(2 * b/(a + b))
                !!Omega1(i, j) = 2 * sqrt(a + b) * (elliptic_inc_ek(par_pi/2, r) - elliptic_inc_ek(0.0_8, r))
                

            end do
        end do
        !$omp end do

		!$omp end parallel

    end subroutine Tab1_Set

end module Tab1


