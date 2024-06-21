module MK
    ! В этом модуде будет храниться функции Монте-Карло, 
    ! а также сечения HH - столкновений
	USE STORAGE
    USE OMP_LIB
    implicit none
    integer, PARAMETER :: MK_n_phi = 768      !  Количество разбиений сечения по углу
    integer, PARAMETER :: MK_n_g = 1301       !  Количество разбиений сечения по относительной скорости
    real(8), PARAMETER :: MK_norm_g = 0.110093       !!  Нормировка скорости (километр в секунду в безразмерном виде)
    real(8), PARAMETER :: MK_norm_sig = 1.0       !!  Нормировка сечения (параметр самой физической задачи)
    integer, PARAMETER :: MK_n_FG = 2000       !!  Точность для розыгрыша g (обратная к Fg)
    real(8) :: MK_norm_A = 0.0       !!  Нормировка функции распределения для вычисления интеграла
    real(8), allocatable :: MK_phi(:)  ! MK_n_phi
    real(8), allocatable :: MK_g(:)   ! MK_n_g
    real(8), allocatable :: MK_sigHH(:, :) ! (MK_n_phi, MK_n_g)
    real(8), allocatable :: MK_sigHH_g(:) ! (MK_n_g)   ! Сечение, проинтегрированное по углу
    real(8), allocatable :: MK_obr_FG(:) ! (MK_n_FG)   ! Обратная функция чтобы найти g от xi
    real(8), allocatable :: MK_obr_Fchi(:, :) ! (MK_n_FG, MK_n_g)   ! Обратная функция чтобы найти g от xi

    real(8), allocatable :: Omega2(:, :, :, :)  !! Сначала две для VH потом две для Vp

    integer(4) :: sensor(3, 2, MK_n_potok)  !(3, 2, :  число потоков)  ! датчики случайных чисел 
        ! Каждому потоку по два датчика


    contains

    function MK_Get_sigma_g(g)
        ! Получить значение сечения в этой точке (без интерполяции)
        real(8), intent(in) :: g
        integer :: i, ng
        real(8) :: MK_Get_sigma_g

        ng = MK_n_g
        do i = 1, MK_n_g
            if(g < MK_g(i)) then
                ng = i
                EXIT
            end if
        end do

        MK_Get_sigma_g = MK_sigHH_g(ng)

    end function MK_Get_sigma_g

    subroutine MK_Read_Sig()
        logical :: exists
        integer(4) :: i, a, b, c, j, k, n1, kk, n3, ii, jj
        real(8) :: x1, x2, x3, S, dchi, g, dg, SS, chi, dphi, phi, u, r, Vx, Vr, Wx, Wr
        character(len=1) :: name


        print*, "MK_Read_Sig -> START" 
        
        allocate(MK_phi(MK_n_phi))
        allocate(MK_g(MK_n_g))
        allocate(MK_sigHH(MK_n_phi, MK_n_g))
        allocate(MK_sigHH_g(MK_n_g))
        allocate(MK_obr_FG(MK_n_FG))
        allocate(MK_obr_Fchi(MK_n_FG, MK_n_g))
        ALLOCATE(Omega2(f1%par_nv1, f1%par_nv2, f1%par_nv1, f1%par_nv2))

		Omega2 = 0.0
        MK_obr_FG = 0.0
        MK_obr_Fchi = par_pi

        if(.True.) then !! Считываем датчики случайных чисел
            inquire(file="rnd_my.txt", exist=exists)
                
            if (exists == .False.) then
                pause "net faila!!!  345434wertew21313edftr3e"
                STOP "net faila!!!"
            end if

            open(1, file = "rnd_my.txt", status = 'old')
                
            do i = 1, MK_n_potok + 5
                read(1,*) a, b, c
                sensor(:, 1, i) = (/ a, b, c /)
                read(1,*) a, b, c
                sensor(:, 2, i) = (/ a, b, c /)
            end do
            
            close(1)
        end if

        if(.True.) then !! Считываем дифференциальное сечение

            n1 = 1

            do j = 1, 8
                write(unit=name,fmt='(i1.1)') j
                inquire(file="all_" // name //".txt", exist=exists)
                    
                if (exists == .False.) then
                    pause "net faila all_1.txt!  gdrgefere"
                    STOP "net faila!!!"
                end if

                open(1, file = "all_" // name //".txt", status = 'old')

                read(1,*) a, b

                do i = 1, a
                    read(1,*) x1
                    MK_g(n1) = x1 * 2.0 * MK_norm_g  !! Здесь надо умножить на 2, так как сечение дано для относительной скорости центра масс (что есть половина от обычной)
                    do k = 1, b
                        read(1,*) x2, x3
                        MK_phi(k) = x2
                        MK_sigHH(k, n1) = x3 * MK_norm_sig/(2.0 * par_pi)
                    end do
                    n1 = n1 + 1
                end do

                close(1)
            end do
        end if

        if(.True.) then  !! Интегрируем сечение по углу
            MK_sigHH_g = 0.0
            do i = 1, MK_n_g
                do j = 1, MK_n_phi
                    if(j > 1 .and. j < MK_n_phi) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0 - (MK_phi(j) + MK_phi(j - 1))/2.0
                    else if( j == 1) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0
                    else
                        dchi = par_pi - (MK_phi(j) + MK_phi(j - 1))/2.0
                    end if
                    MK_sigHH_g(i) = MK_sigHH_g(i) + dchi * MK_sigHH(j, i)
                end do
            end do

            open(1, file = 'sigma_g.txt')

            do i = 1, MK_n_g
                write(1,*) MK_g(i), MK_sigHH_g(i)
            end do

            close(1)

        end if

        if(.True.) then !! Вычисляем A (нормировка функции распределения)
            S = 0.0
            do i = 1, MK_n_g
                g = MK_g(i)
                
                if(g > par_g_max) EXIT

                if(i > 1 .and. i < MK_n_phi) then
                    dg = (MK_g(i + 1) + MK_g(i))/2.0 - (MK_g(i) + MK_g(i - 1))/2.0
                else if( i == 1) then
                    dg = (MK_g(i + 1) + MK_g(i))/2.0
                else
                    dg = MK_g(i) - MK_g(i - 1)
                end if

                do j = 1, MK_n_phi
                    if(j > 1 .and. j < MK_n_phi) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0 - (MK_phi(j) + MK_phi(j - 1))/2.0
                    else if( j == 1) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0
                    else
                        dchi = par_pi - (MK_phi(j) + MK_phi(j - 1))/2.0
                    end if

                    S = S + g**3 * MK_sigHH(j, i) * dchi * dg
                end do
            end do

            MK_norm_A = S * 8.0 * par_pi**2 
            print*, "MK_norm_A = ", MK_norm_A
        end if

        if(.True.) then !! Вычисляем первообразную Fg
            S = 0.0
            do i = 1, MK_n_g
                g = MK_g(i)

                if(g > par_g_max) EXIT

                if(i > 1 .and. i < MK_n_phi) then
                    dg = (MK_g(i + 1) + MK_g(i))/2.0 - (MK_g(i) + MK_g(i - 1))/2.0
                else if( i == 1) then
                    dg = (MK_g(i + 1) + MK_g(i))/2.0
                else
                    dg = MK_g(i) - MK_g(i - 1)
                end if

                S = S + g**3 * MK_sigHH_g(j) * dg
            end do

            MK_obr_FG = par_g_max
            k = 1
            SS = 0.0
            do i = 1, MK_n_g
                g = MK_g(i)

                if(g > par_g_max) EXIT

                if(i > 1 .and. i < MK_n_phi) then
                    dg = (MK_g(i + 1) + MK_g(i))/2.0 - (MK_g(i) + MK_g(i - 1))/2.0
                else if( i == 1) then
                    dg = (MK_g(i + 1) + MK_g(i))/2.0
                else
                    dg = MK_g(i) - MK_g(i - 1)
                end if

                SS = SS + g**3 * MK_sigHH_g(j) * dg

                do while(k * 1.0/MK_n_FG < SS/S)
                    MK_obr_FG(k) = g
                    k = k + 1
                end do
            end do

            open(1, file = 'obr_F_g.txt')

            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_FG(i)
            end do

            close(1)

        end if

        if(.True.) then  !! Вычисляем первообразные для chi при каждом g
            do k = 1, MK_n_g
                g = MK_g(k)
                if(g > par_g_max) EXIT
                S = 0.0
                do j = 1, MK_n_phi
                    if(j > 1 .and. j < MK_n_phi) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0 - (MK_phi(j) + MK_phi(j - 1))/2.0
                    else if( j == 1) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0
                    else
                        dchi = par_pi - (MK_phi(j) + MK_phi(j - 1))/2.0
                    end if

                    S = S + MK_sigHH(j, k) * dchi
                end do
                SS = 0.0

                kk = 1
                do j = 1, MK_n_phi
                    chi = MK_phi(j)
                    if(j > 1 .and. j < MK_n_phi) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0 - (MK_phi(j) + MK_phi(j - 1))/2.0
                    else if( j == 1) then
                        dchi = (MK_phi(j + 1) + MK_phi(j))/2.0
                    else
                        dchi = par_pi - (MK_phi(j) + MK_phi(j - 1))/2.0
                    end if

                    SS = SS + MK_sigHH(j, k) * dchi

                    do while(kk * 1.0/MK_n_FG < SS/S)
                        MK_obr_Fchi(kk, k) = chi
                        kk = kk + 1
                    end do
                end do
            end do

            open(1, file = 'obr_F_chi_100.txt')
            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_Fchi(i, 100)
            end do
            close(1)

            open(1, file = 'Sig_chi_100.txt')
            do i = 1, MK_n_phi
                write(1,*) MK_phi(i), MK_sigHH(i, 100)
            end do
            close(1)

            open(1, file = 'obr_F_chi_200.txt')
            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_Fchi(i, 200)
            end do
            close(1)

            open(1, file = 'Sig_chi_200.txt')
            do i = 1, MK_n_phi
                write(1,*) MK_phi(i), MK_sigHH(i, 200)
            end do
            close(1)

            open(1, file = 'obr_F_chi_300.txt')
            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_Fchi(i, 300)
            end do
            close(1)

            open(1, file = 'Sig_chi_300.txt')
            do i = 1, MK_n_phi
                write(1,*) MK_phi(i), MK_sigHH(i, 300)
            end do
            close(1)

            open(1, file = 'obr_F_chi_400.txt')
            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_Fchi(i, 400)
            end do
            close(1)

            open(1, file = 'Sig_chi_400.txt')
            do i = 1, MK_n_phi
                write(1,*) MK_phi(i), MK_sigHH(i, 400)
            end do
            close(1)

            open(1, file = 'obr_F_chi_500.txt')
            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_Fchi(i, 500)
            end do
            close(1)

            open(1, file = 'Sig_chi_500.txt')
            do i = 1, MK_n_phi
                write(1,*) MK_phi(i), MK_sigHH(i, 500)
            end do
            close(1)

            open(1, file = 'obr_F_chi_600.txt')
            do i = 1, MK_n_FG
                write(1,*) i * 1.0/MK_n_FG, MK_obr_Fchi(i, 600)
            end do
            close(1)

            open(1, file = 'Sig_chi_600.txt')
            do i = 1, MK_n_phi
                write(1,*) MK_phi(i), MK_sigHH(i, 600)
            end do
            close(1)

        end if
        
        if(.False.) then  !! Вычисляем Omega2

            n3 = 200
            dphi = par_pi/n3

            !$omp parallel
            !$omp do private(j, k, phi, S, a, b, u, r, Vx, Vr, Wx, Wr, ii, jj) schedule(dynamic, 1)
            do i = 1, f1%par_nv1
                print*, i, " from ", f1%par_nv1
                call Get_param_Vx(f1, i, Vx)
                do j = 1, f1%par_nv2
                    call Get_param_Vr(f1, j, Vr)
                    do ii = 1, f1%par_nv1
                        call Get_param_Vx(f1, ii, Wx)
                        do jj = 1, f1%par_nv2
                            call Get_param_Vr(f1, jj, Wr)
                            S = 0.0
                            a = (Vx - Wx)**2 + Vr**2 + Wr**2
                            b = 2.0 * Vr * Wr
                            do k = 1, n3
                                phi = (k - 0.5) * dphi
                                u = sqrt(a - b * cos(phi))
                                S = S + u * MK_Get_sigma_g(u) * dphi
                            end do
                            Omega2(i, j, ii, jj) = S * 2.0 * par_pi
                        end do
                    end do
                end do
            end do
            !$omp end do
            !$omp end parallel

            open(1, file = "save_Omega2.bin", FORM = 'BINARY')
            write(1) f1%par_nv1, f1%par_nv2
            write(1) Omega2
            close(1)
        else
            inquire(file= "save_Omega2.bin", exist=exists)
    
            if (exists == .False.) then
                pause "net faila!!!"
                STOP "net faila!!!"
            end if

            open(1, file = "save_Omega2.bin", FORM = 'BINARY', ACTION = "READ")
        
            read(1) ii, jj
            if(ii /= f1%par_nv1 .or. jj /= f1%par_nv2) then
                pause "net faila wfcwvvfrscsdrfescfesrgsr!!!"
                STOP "net faila!!!"
            end if
            read(1) Omega2
            close(1)
        end if

        print*, "MK_Read_Sig -> END" 

    end subroutine MK_Read_Sig

    function MK_g_Get(ksi)
        ! Получить значение сечения в этой точке (без интерполяции)
        real(8), intent(in) :: ksi
        real(8) :: dksi, MK_g_Get
        integer :: n

        dksi = 1.0/MK_n_FG

        n = int(ksi/dksi)
        if(n < 1) n = 1
        if(n > MK_n_FG) n = MK_n_FG

        MK_g_Get = MK_obr_FG(n)

    end function MK_g_Get

    function MK_chi_Get(ksi, g)
        ! Получить значение сечения в этой точке (без интерполяции)
        real(8), intent(in) :: ksi, g
        real(8) :: dksi, MK_chi_Get
        integer :: n, ng, i

        ng = MK_n_g
        do i = 1, MK_n_g
            if(MK_g(i) > g) then
                ng = i
                EXIT
            end if
        end do

        dksi = 1.0/MK_n_FG

        n = int(ksi/dksi)
        if(n < 1) n = 1
        if(n > MK_n_FG) n = MK_n_FG

        MK_chi_Get = MK_obr_Fchi(n, ng)

    end function MK_chi_Get

    function MK_Get(g, chi)
        ! Получить значение сечения в этой точке (без интерполяции)
        real(8), intent(in) :: g, chi
        integer :: i, j, ng, nchi
        real(8) :: MK_Get

        ng = MK_n_g
        do i = 1, MK_n_g
            if(g < MK_g(i)) then
                ng = i
                EXIT
            end if
        end do

        nchi = MK_n_phi
        do i = 1, MK_n_phi
            if(chi < MK_phi(i)) then
                nchi = i
                EXIT
            end if
        end do

        MK_Get = MK_sigHH(nchi, ng)

    end function MK_Get

    subroutine MK_exch_velos(wx, wy, gg, gx, gy, gz, chi, eps, wwx, wwr, wwx1, wwr1)
        ! Считает новые скорости после столкновения (чисто геометрия)
        real(8), intent(in) :: wx, wy, gg, gx, gy, gz, chi, eps
        real(8), intent(out) :: wwx, wwr, wwx1, wwr1
        real(8) :: lx, ly, lz, ksiax, ksiay, ksiaz, gxy, wwy, wwz, wwy1, wwz1

        ksiax = 0.5 * (2.0 * wx - gx)
        ksiay = 0.5 * (2.0 * wy - gy)
        ksiaz = 0.5 * (- gz)

        gxy = sqrt(gx**2 + gy**2)

        lx = gx/gg * cos(chi) - gx * gz/gg/gxy * cos(eps) * sin(chi) + gy/gxy * sin(eps) * sin(chi)
        ly = gy/gg * cos(chi) - gy * gz/gg/gxy * cos(eps) * sin(chi) - gx/gxy * sin(eps) * sin(chi)
        lz = gz/gg * cos(chi) + gxy/gg * cos(eps) * sin(chi)

        wwx = ksiax + 0.5 * gg * lx
        wwy = ksiay + 0.5 * gg * ly
        wwz = ksiaz + 0.5 * gg * lz

        wwx1 = ksiax - 0.5 * gg * lx
        wwy1 = ksiay - 0.5 * gg * ly
        wwz1 = ksiaz - 0.5 * gg * lz

        !print*, "l = ", lx, ly, lz, ksiax, ksiay, ksiaz

        wwr = sqrt(wwy**2 + wwz**2)
        wwr1 = sqrt(wwy1**2 + wwz1**2)
        
        return
    end subroutine MK_exch_velos

    subroutine Test_int()
        integer :: N, i1, i2, i3, i4, i5
        integer :: N1, N2, N3, N4, N5
        real(8) :: S, wx, wr, ksi1, ksi2, ksi3, ksi4, ksi5, eps, phi, chi, g, the
        real(8) :: wwx, wwr, wwx1, wwr1, d1, d2, d3, d4, d5

        wx = 1.0
        wr = 0.2
        S = 0.0
        do N = 1, 1000000
            call M_K_rand(sensor(1, 1, 1), sensor(2, 1, 1), sensor(3, 1, 1), ksi1)
            call M_K_rand(sensor(1, 2, 1), sensor(2, 2, 1), sensor(3, 2, 1), ksi2)
            call M_K_rand(sensor(1, 2, 1), sensor(2, 2, 1), sensor(3, 2, 1), ksi3)
            call M_K_rand(sensor(1, 1, 1), sensor(2, 1, 1), sensor(3, 1, 1), ksi4)
            call M_K_rand(sensor(1, 2, 1), sensor(2, 2, 1), sensor(3, 2, 1), ksi5)

            eps = 2.0 * par_pi * ksi1
            phi = 2.0 * par_pi * ksi2
            the = acos(1.0 - 2.0 * ksi3)

            g = MK_g_Get(ksi4)
            chi = MK_chi_Get(ksi5, g)

            call MK_exch_velos(wx, wr, g, g * sin(the) * cos(phi), g  * sin(the) * sin(phi), g * cos(the), chi, eps, wwx, wwr, wwx1, wwr1)
            ! print*, "______________"
            ! print*, chi, eps
            ! print*, "____"
            ! print*, wwx, wwr, wwx1, wwr1
            ! pause
            S = S + f_maxwell(wwx, wwr, -0.1_8, 0.2_8) * f_maxwell(wwx1, wwr1, -0.1_8, 0.2_8) 
        end do

        S = S * MK_norm_A/1000000

        print*, "1 Test int = ", S

        S = 0.0
        N1 = MK_n_phi
        N2 = 10
        N3 = MK_n_g
        N4 = 10
        N5 = 80

        d2 = 2.0 * par_pi/N2
        d4 = 2.0 * par_pi/N4
        d5 = par_pi/N5

        !!! даже один 5-мерный интеграл численно взять не получается !!!
        ! do i3 = 1, N3
        !     g = MK_g(i3)
        !     if(g > par_g_max) EXIT

        !     if(i3 > 1 .and. i3 < MK_n_phi) then
        !         d3 = (MK_g(i3 + 1) + MK_g(i3))/2.0 - (MK_g(i3) + MK_g(i3 - 1))/2.0
        !     else if( i3 == 1) then
        !         d3 = (MK_g(i3 + 1) + MK_g(i3))/2.0
        !     else
        !         d3 = MK_g(i3) - MK_g(i3 - 1)
        !     end if

        !     do i1 = 1, N1
        !         chi = MK_phi(i1)

        !         if(i1 > 1 .and. i1 < MK_n_phi) then
        !             d1 = (MK_phi(i1 + 1) + MK_phi(i1))/2.0 - (MK_phi(i1) + MK_phi(i1 - 1))/2.0
        !         else if( i1 == 1) then
        !             d1 = (MK_phi(i1 + 1) + MK_phi(i1))/2.0
        !         else
        !             d1 = par_pi - (MK_phi(i1) + MK_phi(i1 - 1))/2.0
        !         end if

        !         do i2 = 1, N2
        !             eps = i2 * d2
        !             do i4 = 1, N4
        !                 phi = i4 * d4
        !                 do i5 = 1, N5
        !                     the = i5 * d5
        !                     call MK_exch_velos(wx, wr, g, g * sin(the) * cos(phi), g  * sin(the) * sin(phi), g * cos(the), chi, eps, wwx, wwr, wwx1, wwr1)
        !                     S = S + g**3 * &
        !                         MK_sigHH(i1, i3) * sin(the) * d1 * d3  
        !                     ! f_maxwell(wwx, wwr, -0.1_8, 0.2_8) * f_maxwell(wwx1, wwr1, -0.1_8, 0.2_8)
        !                 end do
        !             end do
        !         end do
        !     end do
        ! end do

        do i3 = 1, N3
            g = MK_g(i3)
            if(g > par_g_max) EXIT

            if(i3 > 1 .and. i3 < MK_n_phi) then
                d3 = (MK_g(i3 + 1) + MK_g(i3))/2.0 - (MK_g(i3) + MK_g(i3 - 1))/2.0
            else if( i3 == 1) then
                d3 = (MK_g(i3 + 1) + MK_g(i3))/2.0
            else
                d3 = MK_g(i3) - MK_g(i3 - 1)
            end if

            do i1 = 1, N1
                chi = MK_phi(i1)

                if(i1 > 1 .and. i1 < MK_n_phi) then
                    d1 = (MK_phi(i1 + 1) + MK_phi(i1))/2.0 - (MK_phi(i1) + MK_phi(i1 - 1))/2.0
                else if( i1 == 1) then
                    d1 = (MK_phi(i1 + 1) + MK_phi(i1))/2.0
                else
                    d1 = par_pi - (MK_phi(i1) + MK_phi(i1 - 1))/2.0
                end if
                do i5 = 1, N5
                    the = i5 * d5
                    !call MK_exch_velos(wx, wr, g, g * sin(the) * cos(phi), g  * sin(the) * sin(phi), g * cos(the), chi, eps, wwx, wwr, wwx1, wwr1)
                    S = S + g**3 * &
                        MK_sigHH(i1, i3) * sin(the) * d1 * d3  
                    ! f_maxwell(wwx, wwr, -0.1_8, 0.2_8) * f_maxwell(wwx1, wwr1, -0.1_8, 0.2_8)
                end do
            end do
        end do

        S = S * d5 * 2.0 * par_pi * 2.0 * par_pi
        print*, "2 Test int = ", S

    end subroutine Test_int

    subroutine M_K_rand(s1, s2, s3, b)
		integer(4), intent(in out) :: s1, s2, s3
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

end module MK