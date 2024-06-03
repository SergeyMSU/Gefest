
module My_func
    implicit none 
    
    real(8), parameter :: MF_par_pi = acos(-1.0_8) 
    real(8), parameter :: MF_meDmp = (1.0_8/1836.15_8)   ! ќтношение массы электрона к массе протона

    contains

    subroutine get_bazis(ex, ey, ez)
            ! ѕо заданному вектору ex подбирает ey, ez ему перпендикул€рные, образующие правую тройку
            !! ex - должен быть единичным
            real(8), intent(in) :: ex(3)
            real(8), intent(out) :: ey(3)
            real(8), intent(out) :: ez(3)
            real(8) :: norm

            ez(1) = 1.0
            ez(2) = 0.0
            ez(3) = 0.0

            if( 1.0 - dabs(DOT_PRODUCT(ex, ez)) < 0.001) then
                ez(1) = 0.0
                ez(2) = 1.0
                ez(3) = 0.0
            end if

            ey(1) = ex(2) * ez(3) - ex(3) * ez(2)
            ey(2) = ex(3) * ez(1) - ex(1) * ez(3)
            ey(3) = ex(1) * ez(2) - ex(2) * ez(1)
            norm = sqrt(ey(1)**2 + ey(2)**2 + ey(3)**2)
            ey = ey/norm

            ez(1) = ex(2) * ey(3) - ex(3) * ey(2)
            ez(2) = ex(3) * ey(1) - ex(1) * ey(3)
            ez(3) = ex(1) * ey(2) - ex(2) * ey(1)
            norm = sqrt(ez(1)**2 + ez(2)**2 + ez(3)**2)
            ez = ez/norm
        end subroutine get_bazis

    !@cuf attributes(host, device) & 
    integer(4) pure function signum(x)
        implicit none
        real(8), intent(in) :: x
        
        if (x > 0) then
            signum = 1
            return
        else if (x < 0) then
            signum = -1
            return
        else 
            signum = 0
            return
        end if
    end function signum
        
    !@cuf attributes(host, device) & 
    real(8) pure function minmod(x, y)
        implicit none
        real(8), intent(in) :: x, y
        
        if (signum(x) + signum(y) == 0) then
            minmod = 0.0_8
            return
        else
            minmod = ((signum(x) + signum(y)) / 2.0) * min(dabs(x), dabs(y)) ! minmod
            return   
        end if
    end function minmod

    !@cuf attributes(host, device) & 
    real(8) pure function  linear(x1, t1, x2, t2, x3, t3, y)
        ! √лавное значение с параметрами 2
        ! —троим линии между 1 и 2,  2 и 3, потом находим минмодом значение в y
        implicit none
        real(8), intent(in) :: x1, x2, x3, y, t1, t2, t3
        real(8) :: d

        d = minmod((t1 - t2) / (x1 - x2), (t2 - t3) / (x2 - x3))
        linear =  (d * (y - x2) + t2)
        return
    end function linear

    !@cuf attributes(host, device) & 
    real(8) pure function  linear1(x1, t1, x2, t2, y)
        implicit none
        real(8), intent(in) :: x1, x2, y, t1, t2
        real(8) :: d

        d = (t1 - t2) / (x1 - x2)
        linear1 =  (d * (y - x2) + t2)
        return
        
    end function linear1

    !@cuf attributes(host, device) & 
    !real(8) pure function  snos_istoch(r1, t1, r2, t2, r)
    real(8) function  snos_istoch(r1, t1, r2, t2, r)
        ! —нос в источнике пропорчионально степенной функции
        implicit none
        real(8), intent(in) :: r1, t1, r2, t2, r
        real(8) :: al

        al = log(t2/t1) / log(r1/r2)
        snos_istoch = (r1/r)**al * t1
        if(snos_istoch <= 0.0) then
            print*, "ERROR fwvrewvewvrsvgrvr", snos_istoch, al, t1, t2, r1, r2, r
            STOP
        end if
        return
        
    end function snos_istoch



    !@cuf attributes(host, device) & 
    real(8) pure function polar_angle(x, y)
        real(8), intent(in) :: x, y

        if (dabs(x) + dabs(y) < 0.000001) then
            polar_angle = 0.0_8
            return
        end if


        if (x < 0) then
            polar_angle = atan(y / x) + 1.0 * MF_par_pi
        elseif (x > 0 .and. y >= 0) then
            polar_angle = atan(y / x)
        elseif (x > 0 .and. y < 0) then
            polar_angle = atan(y / x) + 2.0 * MF_par_pi
        elseif (y > 0 .and. x >= 0 .and. x <= 0) then
            polar_angle = MF_par_pi / 2.0
        elseif (y < 0 .and. x >= 0 .and. x <= 0) then
            polar_angle =  3.0 * MF_par_pi / 2.0
        end if

        return
    end function polar_angle

    !@cuf attributes(host, device) & 
    subroutine polyar_skorost(phi, Vy, Vz, Vr, Vphi)
        ! ѕеревод скорости из декартовой в пол€рную с помощью пол€рного угла
        ! Variables
        implicit none
        real(8), intent(in) :: phi, Vy, Vz
        real(8), intent(out) :: Vr, Vphi

        Vr = Vy * cos(phi) + Vz * sin(phi)
        Vphi = Vz * cos(phi) - Vy * sin(phi)
    end subroutine polyar_skorost

    !@cuf attributes(host, device) & 
    subroutine dekard_polyar_skorost(phi, Vr, Vphi, Vy, Vz)
        ! ѕеревод скорости из пол€рной в декартову
        ! Variables
        implicit none
        real(8), intent(in) :: phi, Vr, Vphi 
        real(8), intent(out) :: Vy, Vz

        Vy = Vr * cos(phi) - Vphi * sin(phi)
        Vz = Vr * sin(phi) + Vphi * cos(phi)

    end subroutine dekard_polyar_skorost

    pure function matinv4(A) result(B)
        !! Performs a direct calculation of the inverse of a 4?4 matrix.
        real(8), intent(in) :: A(4,4)   !! Matrix
        real(8)             :: B(4,4)   !! Inverse matrix
        real(8)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = &
        1/(A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
        - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
        + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
        - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))

        ! Calculate the inverse of the matrix
        B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
        B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
        B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
        B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
        B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
        B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
        B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
        B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
        B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
        B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
        B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
        B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
        B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
        B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
    end function

    pure function matinv3(A) result(B)
        !! Performs a direct calculation of the inverse of a 3?3 matrix.
        real(8), intent(in) :: A(3,3)   !! Matrix
        real(8)             :: B(3,3)   !! Inverse matrix
        real(8)             :: detinv

        ! Calculate the inverse determinant of the matrix
        detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

        ! Calculate the inverse of the matrix
        B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
        B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
        B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
        B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
        B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
        B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
        B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
        B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
        B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
    end function

    subroutine Print_matrix_real(A)
        real(8), intent(in) :: A(:, :)   !! Matrix
        integer(4) :: mi, ni

        do mi = 1, size(A(:, 1))
            do ni = 1, size(A(1, :))
                write(*,"(F8.2,$)") A(mi,ni)
            end do
            write (*,*) ''
        end do
    end subroutine Print_matrix_real

    !@cuf attributes(host, device) & 
    subroutine dekard_skorost(z, x, y, Vr, Vphi, Vtheta, Vz, Vx, Vy)
        implicit none
        real(8), intent(in) :: x, y, z,  Vr, Vphi, Vtheta
        real(8), intent(out) :: Vx, Vy, Vz
        real(8) :: r_2, the_2, phi_2

        r_2 = sqrt(x * x + y * y + z * z);
        the_2 = acos(z / r_2);
        phi_2 = polar_angle(x, y);
        
        !print*, r_2, the_2, phi_2, Vr, Vphi, Vtheta
        !print*, sin(the_2), cos(phi_2), cos(the_2), sin(phi_2)

        Vx = Vr * sin(the_2) * cos(phi_2) + Vtheta * cos(the_2) * cos(phi_2) - Vphi * sin(phi_2);
        Vy = Vr * sin(the_2) * sin(phi_2) + Vtheta * cos(the_2) * sin(phi_2) + Vphi * cos(phi_2);
        Vz = Vr * cos(the_2) - Vtheta * sin(the_2);
	end subroutine dekard_skorost

    real(8) pure function kvv(x, y, z)
        real(8), intent (in) :: x, y, z
        kvv = x**2 + y**2 + z**2
    end function kvv

    !@cuf attributes(host, device) & 
    subroutine spherical_skorost(z, x, y, Vz, Vx, Vy, Vr, Vphi, Vtheta)
        implicit none
        real(8), intent(in) :: x, y, z, Vx, Vy, Vz
        real(8), intent(out) :: Vr, Vphi, Vtheta
        real(8) :: r_1, the_1, phi_1

        r_1 = sqrt(x * x + y * y + z * z)
        the_1 = acos(z / r_1)
        phi_1 = polar_angle(x, y)

        Vr = Vx * sin(the_1) * cos(phi_1) + Vy * sin(the_1) * sin(phi_1) + Vz * cos(the_1);
        Vtheta = Vx * cos(the_1) * cos(phi_1) + Vy * cos(the_1) * sin(phi_1) - Vz * sin(the_1);
        Vphi = -Vx * sin(phi_1) + Vy * cos(phi_1);
	end subroutine spherical_skorost

    subroutine Sootnosheniya(rho, p, rho_He, rho_Pui, T_Pui, al, rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E)
        ! ‘ункци€, определ€юща€ температуры и концентрации гели€, пикапов и т.д.
        ! al - это зар€д гели€
        ! если al = 1 то вне гелиопаузы
        ! если al = 2, то внутри гелиопаузы
        ! Th - термальные протоны, He - гелий, Pui - пикапы, E - электроны
        ! без параметров, это общие (те, что считаютс€ в ћ√ƒ)
        implicit none
        real(8), intent(in) :: rho, p, rho_He, T_Pui, rho_Pui
        integer(4), intent(in) :: al
        real(8), intent(out), optional :: rho_Th, rho_E, p_Th, p_Pui, T_Th, T_E

        if(PRESENT(rho_Th)) then
            rho_Th = -(MF_meDmp * al * rho_He + 4.0 * (-rho + rho_He))/(4.0 * (1.0 + MF_meDmp)) - rho_Pui
        end if

        if(PRESENT(rho_E)) then
            rho_E = MF_meDmp * (4.0 * rho + (-4.0 + al) * rho_He)/(4.0 * (1.0 + MF_meDmp))
        end if

        if(PRESENT(p_Th)) then
            p_Th = (p - rho_Pui * T_Pui) * (-4.0 * rho + (4.0 + al * MF_meDmp) * rho_He + 4.0 * (1.0 + MF_meDmp) * rho_Pui)/&
            (-8.0 * rho + (7.0 - al + (-1.0 + al) * MF_meDmp) * rho_He + 4.0 * (1.0 + MF_meDmp) * rho_Pui)
        end if

        if(PRESENT(p_Pui)) then
            p_Pui = T_pui * rho_Pui
        end if

        if(PRESENT(T_Th)) then
            T_Th = -4.0 * (1.0 + MF_meDmp) * (p - rho_Pui * T_Pui)/&
            (-8.0 * rho + (7.0 - al + (-1.0 + al)*MF_meDmp) * rho_He + 4.0 * (1.0 - MF_meDmp) * rho_Pui)
        end if

        if(PRESENT(T_E)) then
            if(PRESENT(T_Th)) then
                T_E = T_Th
            else
                T_E = -4.0 * (1.0 + MF_meDmp) * (p - rho_Pui * T_Pui)/&
                (-8.0 * rho + (7.0 - al + (-1.0 + al)*MF_meDmp) * rho_He + 4.0 * (1.0 - MF_meDmp) * rho_Pui)
            end if
        end if


        return
	end subroutine Sootnosheniya

end module My_func