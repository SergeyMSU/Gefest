
module Distfunc
	USE STORAGE
	USE My_func
	USE cgod
	USE Tab1
	implicit none
	
	contains
	
	function f_maxwell(x, y, u, cp)
		USE STORAGE
		implicit none
		real(8), intent(in) :: x, y, u, cp
		real(8) :: f_maxwell
		
		f_maxwell = 1.0/(cp**3 * par_pi**1.5) * exp(-(x - u)**2/cp**2 - (y)**2/cp**2)
	end function f_maxwell
	
	subroutine Get_param_Vx(f, i, Vx)
		TYPE (DistF), intent(in) :: f
		integer, intent(in) :: i
		real(8), intent(out) :: Vx
		
		Vx = f%par_Lv1 + (i - 0.5) * (f%par_Rv1 - f%par_Lv1)/(f%par_nv1)
		return
	end subroutine Get_param_Vx
	
	subroutine Get_param_Vr(f, j, Vr)
		TYPE (DistF), intent(in) :: f
		integer, intent(in) :: j
		real(8), intent(out) :: Vr
		
		Vr = f%par_Lv2 + (j - 0.5) * (f%par_Rv2 - f%par_Lv2)/(f%par_nv2)
		return
	end subroutine Get_param_Vr
	
	subroutine Get_param_x(f, k, x)
		TYPE (DistF), intent(in) :: f
		integer, intent(in) :: k
		real(8), intent(out) :: x
		
		x = f%par_L + (k - 0.5) * (f%par_R - f%par_L)/(f%par_n)
		return
	end subroutine Get_param_x
	
	subroutine Get_param(f, i, j, k, Vx, Vr, x)
		TYPE (DistF), intent(in) :: f
		integer, intent(in) :: i, j, k
		real(8), intent(out) :: Vx, Vr, x
		
		call Get_param_Vx(f, i, Vx)
		call Get_param_Vx(f, j, Vr)
		call Get_param_Vx(f, k, x)
		return
	end subroutine Get_param
	
	subroutine Init_Func(f)
		! Инициализации функции распределения
		TYPE (DistF), intent(in out) :: f
		integer :: i, j, k, n
		real(8) :: Vx, Vr, x, tt
		logical:: exists
		n = 0
		allocate(f%DistF(f%par_nv1, f%par_nv2, f%par_n))
		allocate(f%Q1m(f%par_nv1, f%par_nv2, f%par_n))
		allocate(f%Q1p(f%par_n))
		
		f%Q1m = 0.0
		f%Q1p = 0.0
		f%DistF = 0.0_8
		do k = 1, f%par_n
			call Get_param_x(f, k, x)
			do j = 1,  f%par_nv2
				call Get_param_Vr(f, j, Vr)
				do i = 1, f%par_nv1
					call Get_param_Vx(f, i, Vx)
					if(x <= 0.0) f%DistF(i, j, k) = f_maxwell(Vx, Vr, f%par_Usr, f%par_c)
				end do
			end do
		end do

		if(ALLOCATED(time_step) == .False.) then
			inquire(file="time_step.txt", exist=exists)
        
			if (exists == .False.) then
				pause "net faila !!! f fwf234r435234r3rw3r4"
				STOP "net faila!!!"
			end if

			open(1, file = "time_step.txt", ACTION = "READ", status = 'old')
			read(1,*) n
			allocate(time_step(n))
			allocate(step_algoritm(f%par_nv1, n))
			
			do i = 1, n
				read(1,*) tt
				time_step(i) = tt
			end do

			do i = 1, n
				do j = 1, f%par_nv1
					read(1,*) k
					step_algoritm(j, i) = k
				end do
			end do

			! print*, step_algoritm(:, 5)
			! print*, step_algoritm(:, 6)

		end if

		
	end subroutine Init_Func

	subroutine Init_FuncGD(g)
		TYPE (GD), intent(in out) :: g
		integer(4) :: i


		!allocate(g%X(g%par_n))
		allocate(g%par(3, g%par_n - 1, 2))

		! Заполнения начальных условий
		do i = 1, g%par_n - 1
			g%par(1, i, 1) = 1.0_8
			g%par(1, i, 2) = 1.0_8

			g%par(2, i, 1) = 0.0_8
			g%par(2, i, 2) = 0.0_8

			g%par(3, i, 1) = 0.3_8
			g%par(3, i, 2) = 0.3_8
		end do


	end subroutine Init_FuncGD

	function get_coordinat_yzel(i, L, R, g)
		TYPE (GD), intent(in out) :: g
		integer(4), intent(in) :: i
		real(8), intent(in) :: L, R
		real(8) :: get_coordinat_yzel
		
		get_coordinat_yzel = L + (i - 1) * (R - L)/(g%par_n - 1);

	end function get_coordinat_yzel

	subroutine Start_GD(T, g, now2)
		TYPE (GD), intent(in out) :: g
		real(8), intent(in) :: T
		integer(4), intent(in) :: now2  !! Какой номер меняем
		integer(4) :: i, j, now
		real(8) :: dx, dx_do, newL, newR, w
		real(8) :: par1(size(g%par(:, 1, 1)))
		real(8) :: par2(size(g%par(:, 1, 1)))
		real(8) :: POTOK(size(g%par(:, 1, 1)))
		real(8) :: normal, dsl, dsp, dsc
		real(8) :: qqq1(9), qqq2(9), POTOK2(9), ro, p, u, ro2, p2, u2, pp, c0, Mach
		logical :: contact, lin_Mach
		integer(4) :: kdir, idgod, KOBL

		c0 = sqrt(g%par_ggg * 0.3/1.0)
		contact = .False.
		KOBL = 0
		kdir = 0
		idgod = 0

		now = mod(now2, 2) + 1
		g%time_step = 100000.0

		!! Считаем новое положение левой границы
		!newL = g%par_L - 0.1 * T! + (g%par(2, 1, now) - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now))) * T
		! if(g%par(1, 1, now) > 0.00001) then
		! 	newL = g%par_L  + (g%par(2, 1, now) - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now))) * T
		! else
		! 	newL = g%par_L  + (g%par(2, 1, now)) * T
		! end if

		! qqq1(1) = 0.000001
		! qqq1(5) = 0.000001
		! qqq1(2) = g%par(2, 1, now)
		! qqq1(3) = 0.0_8
		! qqq1(4) = 0.0_8
		! qqq1(6:8) = 0.0_8
		! qqq1(9) = 0.0_8

		! qqq2(1) = g%par(1, 1, now)
		! qqq2(5) = g%par(3, 1, now)
		! qqq2(2) = g%par(2, 1, now)
		! qqq2(3) = 0.0_8
		! qqq2(4) = 0.0_8
		! qqq2(6:8) = 0.0_8
		! qqq2(9) = 0.0_8
		! w = 0.0_8
		! call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
		! 						1.0_8, 0.0_8, 0.0_8, 1.0_8, &
		! 						w, qqq1(1:8), qqq2(1:8), &
		! 						dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
		! 						POTOK2, kontact_ = contact)
		! newL = g%par_L + dsl * T


		!newL = g%par_L  - (c0) * T
		Mach = dabs(g%par(2, 1, now))/sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now))
		if(Mach < 7.0) then
			lin_Mach = .False.
			newL = g%par_L  + (g%par(2, 1, now) - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now))) * T
		else
			!g%start = .False.
			lin_Mach = .True.
			newL = g%par_L + (g%par(2, 1, now) + 7.0 * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now)) ) * T
			!newL = g%par_L + (g%par(2, 1, now)) * T
		end if

		newR = g%par_R

		dx = get_coordinat_yzel(2, newL, newR, g) - get_coordinat_yzel(1, newL, newR, g)
		dx_do = get_coordinat_yzel(2, g%par_L, g%par_R, g) - get_coordinat_yzel(1, g%par_L, g%par_R, g)


		do i = 1, g%par_n - 1
			
			par1 = g%par(:, i, now)
			POTOK = 0.0

			do j = 1, 2
				contact = .False.
				dsl = 0.0
				dsp = 0.0
				if(j == 1) then ! Правая грань
					if(i == g%par_n - 1) then
						par2 = par1
					else
						par2 = g%par(:, i + 1, now)
					end if
					normal = 1.0_8
					w = normal * (get_coordinat_yzel(i + 1, newL, newR, g) - get_coordinat_yzel(i + 1, g%par_L, g%par_R, g))/T
				else
					w = (get_coordinat_yzel(i, newL, newR, g) - get_coordinat_yzel(i, g%par_L, g%par_R, g))/T
					if(i == 1) then !! Граница с вакуумом
						if(lin_Mach == .False.) then
							par2(1) = par1(1)
							par2(3) = par1(3)
							par2(2) = -par1(2) + 2.0 * w
						else
							par2(1) = par1(1)
							par2(3) = par1(3)
							par2(2) = par1(2)
						end if
						! par2(1) = 0.000001
						! par2(3) = 0.000001
						! par2(2) = par1(2)
						!print*, "w = ", w
					else
						par2 = g%par(:, i - 1, now)
					end if
					normal = -1.0_8
					w = normal * (get_coordinat_yzel(i, newL, newR, g) - get_coordinat_yzel(i, g%par_L, g%par_R, g))/T
					
				end if

				qqq1(1) = par1(1)
				qqq1(5) = par1(3)
				qqq1(2) = par1(2)
				qqq1(3) = 0.0_8
				qqq1(4) = 0.0_8
				qqq1(6:8) = 0.0_8
				qqq1(9) = 0.0_8

				qqq2(1) = par2(1)
				qqq2(5) = par2(3)
				qqq2(2) = par2(2)
				qqq2(3) = 0.0_8
				qqq2(4) = 0.0_8
				qqq2(6:8) = 0.0_8
				qqq2(9) = 0.0_8
				
				if(j == 2 .and. i == 1 .and. lin_Mach == .False.) contact = .True.

				!if(.False.) then!
				if(j == 2 .and. i == 1 .and. lin_Mach == .False.) then
					POTOK2 = 0.0
				else
					if(j == 2 .and. i == 1 .and. lin_Mach == .True.) w = 0.0_8
					call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
								normal, 0.0_8, 0.0_8, 1.0_8, &
								w, qqq1(1:8), qqq2(1:8), &
								dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
								POTOK2, kontact_ = contact)
				end if

				!$omp critical
					g%time_step = min( g%time_step,  0.2 * dx/( max(dabs(dsl), dabs(dsp)) + dabs(w) ))
				!$omp end critical

				if (idgod == 2) print*, "ERROR 4i0u43h9h43t3r434r3"

				POTOK(1) = POTOK(1) + POTOK2(1)
				POTOK(2) = POTOK(2) + POTOK2(2)
				POTOK(3) = POTOK(3) + POTOK2(5)
			end do

			ro = par1(1)
			p = par1(3)
			u = par1(2)

			ro2 = par1(1) * dx_do/dx - T * (POTOK(1) / dx)
			if(ro2 < 0.0) then
				print*, "ro2 < 0", ro2, ro, i
			end if
			u2 = (ro * u * dx_do/dx - T * ( POTOK(2) / dx )) / ro2
			p2 = ((  ( p / (g%par_ggg - 1.0) + 0.5 * ro * u**2)  * dx_do/dx   &
                        - T * (POTOK(3)/ dx) ) - 0.5 * ro2 * u2**2 ) * (g%par_ggg - 1.0)
			if(p2 < 0.0) then
				print*, "p2 < 0", p2, p, i
				p2 = 0.000001
			end if
			
			g%par(1, i, now2) = ro2
			g%par(3, i, now2) = p2
			g%par(2, i, now2) = u2

		end do

		g%par_L = newL

	end subroutine Start_GD

	subroutine Get_GD(x, g, par, now)
		TYPE (GD), intent(in) :: g
		real(8), intent(in) :: x
		integer(4), intent(in) :: now
		real(8), intent(out) :: par(3)
		integer(4) :: i1, i2
		real(8) :: x1, x2, al

		if(x < g%par_L .or. x > g%par_R) then
			par = 0.0_8
			return
		end if

		i1 = INT((g%par_n - 1) * (x - g%par_L)/(g%par_R - g%par_L)  + 0.5)
		i2 = i1 + 1

        if(i1 < 1) then
            i1 = 1
            i2 = 2
        end if

        if(i2 > g%par_n - 1) then
            i2 = g%par_n - 1
            i1 = i2 - 1
        end if

        if(x > g%par_n) then
            i2 = g%par_n - 1
            i1 = i2 - 1
        end if
		
        x1 = g%par_L + (g%par_R - g%par_L) * (i1 - 0.5)/(g%par_n - 1)
        x2 = g%par_L + (g%par_R - g%par_L) * (i2 - 0.5)/(g%par_n - 1)
		al = (x - x1)/(x2 - x1)
		par = (1.0 - al) * g%par(:, i1, now) + al * g%par(:, i2, now)
	end subroutine Get_GD
	
	subroutine Get_Func(f, i, j, x, ff)
		TYPE (DistF), intent(in) :: f
		integer, intent(in) :: i, j
		real(8), intent(out) :: ff      ! Значение функции в этой точке
		real(8), intent(in) :: x
		integer :: k1, k2
		real(8) :: dd, al, Vx, Vr, L, R
		
		if(x >= f%par_R) then
			ff = 0.0_8
			! call Get_param(f, i, j, 3, Vx, Vr, dd)
			! ff = f_maxwell(Vx, Vr, f%par_Usr, f%par_c)
			return
		end if
		
		if(x <= f%par_L) then
			call Get_param_Vr(f, j, Vr)
			call Get_param_Vx(f, i, Vx)
			ff = f_maxwell(Vx, Vr, f%par_Usr, f%par_c)
			return
		end if
		
		dd = (f%par_R - f%par_L)/f%par_n
		k1 = INT((x - f%par_L)/dd) + 1
		k2 = k1 + 1
		
		L = f%par_L + k1 * dd - dd/2.0
		R = f%par_L + k2 * dd - dd/2.0
		if(L > x) then
			k1 = k1 - 1
			k2 = k2 - 1
			L = f%par_L + k1 * dd - dd/2.0
			R = f%par_L + k2 * dd - dd/2.0
		end if
		al = (x - L)/(R - L)
		
		if(k2 > f%par_n) then
			ff = (1.0_8 - al) * f%DistF(i, j, k1)
			return
		end if
		
		if(k1 < 1) then
			call Get_param_Vr(f, j, Vr)
			call Get_param_Vx(f, i, Vx)
			!ff = (1.0_8 - al) * f_maxwell(Vx, Vr, f%par_Usr, f%par_c) + al * f%DistF(i, j, k2)
			ff = f_maxwell(Vx, Vr, f%par_Usr, f%par_c)
			return
		end if
		
		ff = (1.0_8 - al) * f%DistF(i, j, k1) + al * f%DistF(i, j, k2)
		return
	end subroutine Get_Func

	subroutine Start(T)
		real(8), intent(in out) :: T
		real(8) :: dt, tt, VV, dx
		integer :: step
		tt = 0.0
		dx = (f1%par_R - f1%par_L)/f1%par_n
		!dt = 0.5 * (f1%par_R - f1%par_L)/f1%par_n/f1%par_Rv1
		call Get_param_Vx(f1, f1%par_nv1/2 + 1, VV)
		dt = 1.0_8 * dx/VV
		step = 0
		print*, "Start = dt = ", dt, " - dx VV  ", dx, VV
		DO WHILE (tt < T) 
			step = step + 1
            !call Integrate(f1, f2, dt)
            !call Integrate_Yorming(f1, f2, dt)

			!call Calc_Q(f1, gd1, 1)
            call Integrate_Protiv_potoka2(f1, f2, dt, 1)
			
			!call Integrate(f2, f1, dt)
			!call Integrate_Yorming(f2, f1, dt)

			!call Calc_Q(f2, gd1, 1)
			call Integrate_Protiv_potoka2(f2, f1, dt, 1)
			tt = tt + 2.0 * dt
			if(mod(step, 1) == 0) then
				print*, step, tt, T
			end if
      	END DO 

		T = tt

	end subroutine Start

	subroutine Integrate(ff1, ff2, dt)
		TYPE (DistF), intent(in out) :: ff1, ff2
		real(8), intent(in) :: dt
		integer :: i, j, k
		real(8) :: x0, x1, Vx, Vr, ff, x2, fff


		do k = 1, ff1%par_n
			call Get_param_x(ff1, k, x0)
			do i = 1, ff1%par_nv1
				call Get_param_Vx(ff1, i, Vx)
				x1 = x0 - dt * Vx
				!x2 = x0 - 2.0 * dt * Vx
				do j = 1, ff1%par_nv2
					call Get_Func(ff1, i, j, x1, ff)
					!call Get_Func(ff1, i, j, x2, fff)
					ff2%DistF(i, j, k) = ff
					!ff2%DistF(i, j, k) = 4.0 * ff/3.0 - fff/3.0
				end do
			end do
		end do

	end subroutine Integrate

	subroutine Integrate_Yorming(ff1, ff2, dt)
		TYPE (DistF), intent(in out) :: ff1, ff2
		real(8), intent(in) :: dt
		integer :: i, j, k
		real(8) :: nu, omega, u2p, u2m, u1m, u1, u1p, u1mm, Vx, dx

		dx = (ff1%par_R - ff1%par_L)/ff1%par_n
		!nu = 1.0
		!omega = (4.0_8 * nu**2 + 1.0_8) * (4.0_8 - nu**2)/5.0_8
		do k = 3, ff1%par_n-3
			do i = 1, ff1%par_nv1
				call Get_param_Vx(ff1, i, Vx)
				nu = Vx * dt/dx
				if(dabs(nu) > 1.0) print*, "Error 1"
				!omega = 4.0 * nu**2 - nu**4
				!omega = (4.0 * nu**2 - nu**4 + 3.0)/2.0
				omega = (4.0 * nu**2 + 1.0) * (4.0 - nu**2)/5.0   ! убираем дисперсию
				!print*, omega, " - ", nu, 4.0 * nu**2 - nu**4
				if(omega > 3.0 ) print*, "Error 2"
				if(omega <  4.0 * nu**2 - nu**4) print*, "Error 3"
				!pause
				do j = 1, ff1%par_nv2
					!if(Vx > 0.0) then
						u1m = ff1%DistF(i, j, k - 1) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k) - ff1%DistF(i, j, k - 1))
						u1mm = ff1%DistF(i, j, k - 2) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k - 1) - ff1%DistF(i, j, k - 2))
						u1 = ff1%DistF(i, j, k) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k + 1) - ff1%DistF(i, j, k))
						u1p = ff1%DistF(i, j, k + 1) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k + 2) - ff1%DistF(i, j, k + 1))
						u2m = 0.5_8 * (ff1%DistF(i, j, k - 1) + u1m - 2.0_8/3.0_8 * nu * (u1m - u1mm))
						u2p = 0.5_8 * (ff1%DistF(i, j, k + 1) + u1p - 2.0_8/3.0_8 * nu * (u1p - u1))
						ff2%DistF(i, j, k) = ff1%DistF(i, j, k) - nu/24.0_8 * (-2 * ff1%DistF(i, j, k + 2) + 7 * ff1%DistF(i, j, k + 1) - &
						7 * ff1%DistF(i, j, k - 1) + 2 * ff1%DistF(i, j, k - 2)) - 3.0_8/8.0_8 * nu * (u2p - u2m) - &
						omega/24.0_8 * (ff1%DistF(i, j, k + 2) - 4 * ff1%DistF(i, j, k + 1) + 6 * ff1%DistF(i, j, k) - 4 * &
						ff1%DistF(i, j, k - 1) + ff1%DistF(i, j, k - 2))
						if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
					! else 
					! 	u1m = ff1%DistF(i, j, k - 1) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k) - ff1%DistF(i, j, k - 1))
					! 	u1mm = ff1%DistF(i, j, k - 2) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k - 1) - ff1%DistF(i, j, k - 2))
					! 	u1 = ff1%DistF(i, j, k) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k + 1) - ff1%DistF(i, j, k))
					! 	u1p = ff1%DistF(i, j, k + 1) - 2.0_8/3.0_8 * nu * (ff1%DistF(i, j, k + 2) - ff1%DistF(i, j, k + 1))
					! 	u2m = 0.5_8 * (ff1%DistF(i, j, k - 1) + u1m - 2.0_8/3.0_8 * nu * (u1m - u1mm))
					! 	u2p = 0.5_8 * (ff1%DistF(i, j, k + 1) + u1p - 2.0_8/3.0_8 * nu * (u1p - u1))
					! 	ff2%DistF(i, j, k) = ff1%DistF(i, j, k) - nu/24.0_8 * (-2 * ff1%DistF(i, j, k + 2) + 7 * ff1%DistF(i, j, k + 1) - &
					! 	7 * ff1%DistF(i, j, k - 1) + 2 * ff1%DistF(i, j, k - 2)) - 3.0_8/8.0_8 * nu * (u2p - u2m) - &
					! 	omega/24.0_8 * (ff1%DistF(i, j, k + 2) - 4 * ff1%DistF(i, j, k + 1) + 6 * ff1%DistF(i, j, k) - 4 * &
					! 	ff1%DistF(i, j, k - 1) + ff1%DistF(i, j, k - 2))
					! end if
				end do
			end do
		end do

	end subroutine Integrate_Yorming
	
	subroutine Integrate_Protiv_potoka(ff1, ff2, ff3, TT)
		TYPE (DistF), intent(in out) :: ff1, ff2, ff3
		real(8), intent(in out) :: TT  ! Общее время решения
		integer :: now
		integer :: i, j, k, st
		integer :: nn(ff1%par_nv1)
		real(8) :: dt, Vx_max, dx, nu, Vx, Vr, Time, Vx_min

		st = 1
		do i = ff1%par_nv1/2 + 1, ff1%par_nv1
			nn(i) = st
			st = st + 2
		end do
		st = 1
		do i = ff1%par_nv1/2, 1, -1
			nn(i) = st
			st = st + 2
		end do

		dx = (ff1%par_R - ff1%par_L)/ff1%par_n
		call Get_param_Vx(ff1, ff1%par_nv1, Vx_max)
		call Get_param_Vx(ff1, ff1%par_nv1/2 + 1, Vx_min)
		dt = dx/Vx_min                                    ! Шаг по времени (для функции распределения)
		Time = 0.0


		do while(Time < TT)

			do i = 1, ff1%par_nv1
				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						ff3%DistF(i, j, k) = ff1%DistF(i, j, k)
					end do
				end do
			end do

			do st = 1, ff1%par_nv1 - 1   ! Делаем столько шагов по времени (разбиваем шаги по времени на партии)
				Time = Time + dt
				do i = 1, ff1%par_nv1
					call Get_param_Vx(ff1, i, Vx)
					nu = Vx * dt/dx
					do k = 3, ff1%par_n-3
						do j = 1, ff1%par_nv2
							call Get_param_Vr(ff1, j, Vr)
							
							if(mod(st, ff1%par_nv1 - nn(i)) == 0) then
								if(Vx > 0.0) then
									nu = 1.0
									ff2%DistF(i, j, k) = ff3%DistF(i, j, k) - nu * (ff3%DistF(i, j, k) - ff3%DistF(i, j, k-1)) + &
										0.5_8 * nu * (nu - 1.0_8) * (ff3%DistF(i, j, k) - 2.0_8 * ff3%DistF(i, j, k - 1) + ff3%DistF(i, j, k-2)) !&
										!+ time * proton(1, k) * (-ff1%Q1p(k) * ff1%DistF(i, j, k) + &
										!					f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
								else
									nu = -1.0
									ff2%DistF(i, j, k) = ff3%DistF(i, j, k) + nu * (3.0_8/2.0_8 * ff3%DistF(i, j, k) - &
										2.0_8 * ff3%DistF(i, j, k + 1) + 0.5_8 * ff3%DistF(i, j, k + 2)) + &
										0.5_8 * nu**2 * (ff3%DistF(i, j, k) - 2.0_8 * ff3%DistF(i, j, k + 1) + ff3%DistF(i, j, k+2))! &
										!+ time * proton(1, k) * (-ff1%Q1p(k) * ff1%DistF(i, j, k) + &
										!						f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
								end if
							else
								if(Vx > 0.0) then
									ff2%DistF(i, j, k) = ff1%DistF(i, j, k) - nu * (ff1%DistF(i, j, k) - ff1%DistF(i, j, k-1)) + &
										0.5_8 * nu * (nu - 1.0_8) * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k - 1) + ff1%DistF(i, j, k-2)) !&
										!+ time * proton(1, k) * (-ff1%Q1p(k) * ff1%DistF(i, j, k) + &
										!					f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
								else
									ff2%DistF(i, j, k) = ff1%DistF(i, j, k) + nu * (3.0_8/2.0_8 * ff1%DistF(i, j, k) - &
										2.0_8 * ff1%DistF(i, j, k + 1) + 0.5_8 * ff1%DistF(i, j, k + 2)) + &
										0.5_8 * nu**2 * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k + 1) + ff1%DistF(i, j, k+2))! &
										!+ time * proton(1, k) * (-ff1%Q1p(k) * ff1%DistF(i, j, k) + &
										!						f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
								end if
							end if
							if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
						end do ! j
					end do! k

					if(mod(st, ff1%par_nv1 - nn(i)) == 0) then
						do k = 3, ff1%par_n-3
							do j = 1, ff1%par_nv2
								ff3%DistF(i, j, k) = ff2%DistF(i, j, k)
							end do
						end do
					end if

				end do !  i

				do i = 1, ff1%par_nv1
					do k = 3, ff1%par_n-3
						do j = 1, ff1%par_nv2
							ff1%DistF(i, j, k) = ff2%DistF(i, j, k)
						end do
					end do
				end do

			end do ! st
		end do ! while

		TT = time
	end subroutine Integrate_Protiv_potoka

	subroutine Integrate_Protiv_potoka2(ff1, ff2, dt, now)
		TYPE (DistF), intent(in out), TARGET :: ff1, ff2
		real(8), intent(in) :: dt  ! Здесь dt фиктивное, так как nu фиксировано
		integer, intent(in) :: now  ! Здесь dt фиктивное, так как nu фиксировано
		integer :: i, j, k, st
		real(8) :: nu, omega, u2p, u2m, u1m, u1, u1p, u1mm, Vx, dx, Vr, x
		real(8) :: lf(ff1%par_nv2, 2), time

		TYPE (DistF), POINTER :: ptr_f1, ptr_f2, ptr_swap
		real(8) :: proton(3, ff1%par_n)  !(rho, u, cp)
		real(8) :: par(3)

		!! Заполним параметры протонов (концентрацию)
		do k = 1, ff1%par_n	
			call Get_param_x(ff1, k, x)
			call Get_GD(x, gd1, par, now)
			proton(1, k) = par(1)
			proton(2, k) = par(2)
			if(par(1) > 0.000001) then
				proton(3, k) = sqrt(par(3)/par(1))
			else
				proton(3, k) = 1.0
			end if
		end do

		!ptr_f1 => ff1
		!ptr_f2 => ff2

		dx = (ff1%par_R - ff1%par_L)/ff1%par_n

		!nu = 1.0
		!omega = (4.0_8 * nu**2 + 1.0_8) * (4.0_8 - nu**2)/5.0_8
		nu = 1.0_8

		do st = 1, ff1%par_nv1 - 1
			!call Calc_Q(ptr_f1, gd1, 1)

			nu = 1.0_8
			do i = ff1%par_nv1/2 + 1, ff1%par_nv1  !! Для правой части функции распределения
				
				if(st <= (ff1%par_nv1 - i) * 2) CYCLE

				! print*, i, st
				! pause

				call Get_param_Vx(ff1, i, Vx)
				time = 0.0! dabs(dx/Vx)

				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						call Get_param_Vr(ff1, i, Vr)
						ff2%DistF(i, j, k) = ff1%DistF(i, j, k-1)
						!ff2%DistF(i, j, k) = ff1%DistF(i, j, k) - nu * (ff1%DistF(i, j, k) - ff1%DistF(i, j, k-1)) + &
						!	0.5_8 * nu * (nu - 1.0_8) * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k - 1) + ff1%DistF(i, j, k-2)) !&
							!+ time * proton(1, k) * (-ff1%Q1p(k) * ff1%DistF(i, j, k) + &
							!					f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
					end do
				end do

				!par_ch1(i) = par_ch1(i) + 1
				!if(i==26) print*, "=", ptr_f2%DistF(26, 5, 330)

				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						ff1%DistF(i, j, k) = ff2%DistF(i, j, k)
					end do
				end do
			end do
			! ptr_swap => ptr_f2
			! ptr_f2 => ptr_f1
			! ptr_f1 => ptr_swap

			nu = -1.0_8
			do i = ff1%par_nv1/2, 1, -1  !! Для Левой части функции распределения

				if(st <= (i - 1) * 2) CYCLE

				call Get_param_Vx(ff1, i, Vx)
				time = 0.0! dabs(dx/Vx)

				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						call Get_param_Vr(ff1, i, Vr)
						ff2%DistF(i, j, k) = ff1%DistF(i, j, k) + nu * (3.0_8/2.0_8 * ff1%DistF(i, j, k) - &
						2.0_8 * ff1%DistF(i, j, k + 1) + 0.5_8 * ff1%DistF(i, j, k + 2)) + &
						0.5_8 * nu**2 * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k + 1) + ff1%DistF(i, j, k+2))! &
						!+ time * proton(1, k) * (-ff1%Q1p(k) * ff1%DistF(i, j, k) + &
						!						f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
					end do
				end do

				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						ff1%DistF(i, j, k) = ff2%DistF(i, j, k)
					end do
				end do
				
			end do

		end do



		do i = 1, ff1%par_nv1
			do k = 3, ff1%par_n-3
				do j = 1, ff1%par_nv2
					ff2%DistF(i, j, k) = ff1%DistF(i, j, k)
				end do
			end do
		end do

		! print*, "_____"
		! do j = 1, ff1%par_nv2
		! 	print*, ff2%DistF(26, j, 330)
		! end do

		!print*, ff2%DistF(30, 5, 330), ff1%DistF(30, 5, 330)

	end subroutine Integrate_Protiv_potoka2

	subroutine Integrate_Protiv_potoka3(ff1, ff2, dt, now)
		TYPE (DistF), intent(in out), TARGET :: ff1, ff2
		real(8), intent(in) :: dt  ! Здесь dt фиктивное, так как nu фиксировано
		integer, intent(in) :: now  ! Здесь dt фиктивное, так как nu фиксировано
		integer :: i, j, k, st
		real(8) :: nu, omega, u2p, u2m, u1m, u1, u1p, u1mm, Vx, dx, Vr, x
		real(8) :: lf(ff1%par_nv2, 2), time

		TYPE (DistF), POINTER :: ptr_f1, ptr_f2, ptr_swap
		real(8) :: proton(3, ff1%par_n)  !(rho, u, cp)
		real(8) :: par(3)

		!! Заполним параметры протонов (концентрацию)
		do k = 1, ff1%par_n	
			call Get_param_x(ff1, k, x)
			call Get_GD(x, gd1, par, now)
			proton(1, k) = par(1)
			proton(2, k) = par(2)
			if(par(1) > 0.000001) then
				proton(3, k) = sqrt(par(3)/par(1))
			else
				proton(3, k) = 1.0
			end if
		end do

		ptr_f1 => ff1
		ptr_f2 => ff2

		dx = (ff1%par_R - ff1%par_L)/ff1%par_n

		!nu = 1.0
		!omega = (4.0_8 * nu**2 + 1.0_8) * (4.0_8 - nu**2)/5.0_8
		nu = 1.0_8
		do i = ff1%par_nv1/2 + 1, ff1%par_nv1  !! Для правой части функции распределения
			call Get_param_Vx(ff1, i, Vx)
			time = dx/Vx
			!nu = Vx * dt/dx
			!print*, nu, nu/i, Vx, i, nu/(abs(i - ff1%par_nv1/2) * 2 - 1)
			!pause

			do st = 1, abs(i - ff1%par_nv1/2) * 2 - 1
				!print*, i, st
				!pause
				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						call Get_param_Vr(ff1, i, Vr)
						ptr_f2%DistF(i, j, k) = ptr_f1%DistF(i, j, k) - nu * (ptr_f1%DistF(i, j, k) - ptr_f1%DistF(i, j, k-1)) + &
							0.5_8 * nu * (nu - 1.0_8) * (ptr_f1%DistF(i, j, k) - 2.0_8 * ptr_f1%DistF(i, j, k - 1) + ptr_f1%DistF(i, j, k-2)) !&
							!+ time * proton(1, k) * (-ff1%Q1p(k) * ptr_f1%DistF(i, j, k)! + &
							!					f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
					end do
				end do

				par_ch1(i) = par_ch1(i) + 1
				if(i==26) print*, "=", ptr_f2%DistF(26, 5, 330)
				
				ptr_swap => ptr_f2
				ptr_f2 => ptr_f1
				ptr_f1 => ptr_swap
			end do

			if(mod(abs(i - ff1%par_nv1/2) * 2 - 1, 2) == 0) then
				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						f2%DistF(i, j, k) = f1%DistF(i, j, k)
					end do
				end do
			end if
		end do


		ptr_f1 => ff1
		ptr_f2 => ff2
		nu = -1.0_8

		do i = ff1%par_nv1/2, 1, -1  !! Для Левой части функции распределения
			call Get_param_Vx(ff1, i, Vx)
			time = dx/Vx
			!nu = Vx * dt/dx
			!print*, nu, nu/i, Vx, i, nu/(abs(i - ff1%par_nv1/2) * 2 - 1)
			!pause

			do st = 1, abs(i - ff1%par_nv1/2 - 1) * 2 - 1
				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						call Get_param_Vr(ff1, i, Vr)
						ptr_f2%DistF(i, j, k) = ptr_f1%DistF(i, j, k) + nu * (3.0_8/2.0_8 * ptr_f1%DistF(i, j, k) - &
						2.0_8 * ptr_f1%DistF(i, j, k + 1) + 0.5_8 * ptr_f1%DistF(i, j, k + 2)) + &
						0.5_8 * nu**2 * (ptr_f1%DistF(i, j, k) - 2.0_8 * ptr_f1%DistF(i, j, k + 1) + ptr_f1%DistF(i, j, k+2)) &
						+ time * proton(1, k) * (-ff1%Q1p(k) * ptr_f1%DistF(i, j, k) + &
												f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
					end do
				end do
				
				ptr_swap => ptr_f2
				ptr_f2 => ptr_f1
				ptr_f1 => ptr_swap
			end do

			if(mod(abs(i - ff1%par_nv1/2 - 1) * 2 - 1, 2) == 0) then
				do k = 3, ff1%par_n-3
					do j = 1, ff1%par_nv2
						f2%DistF(i, j, k) = f1%DistF(i, j, k)
					end do
				end do
			end if
		end do

		!print*, "ff = ", ff2%DistF(30, 5, 330), ff1%DistF(30, 5, 330)
		! print*, "_____"
		! do j = 1, ff1%par_nv2
		! 	print*, ff2%DistF(26, j, 330)
		! end do

	end subroutine Integrate_Protiv_potoka3

	subroutine Calc_Q(ff, g, now)
		TYPE (DistF), intent(in out), TARGET :: ff
		TYPE (GD), intent(in) :: g
		integer(4), intent(in) :: now
		integer(4) :: k, i, j, ii, jj
		real(8) :: x, Vx, Vr, Wx, Wr
		real(8) :: dWx, dWr, A, B
		real(8) :: par(3)
		real(8) :: u, S, cp, uz

		dWx = (ff%par_Rv1 - ff%par_Lv1)/(ff%par_nv1)
		dWr = (ff%par_Rv2 - ff%par_Lv2)/(ff%par_nv2)

		do k = 1, ff%par_n
			!print*, k, ff%par_n
			call Get_param_x(ff, k, x)  ! Получаем x координату функции
			call Get_GD(x, g, par, now)
			if(par(1) < 0.0001) then
				ff%Q1m(:, :, k) = 0.0_8
				ff%Q1p(k) = 0.0_8
				CYCLE
			end if

			S = 0.0
			
			! Посчитаем концентрацию водорода, если водорода нет, то и перезаряжать не надо
			do j = 1, ff%par_nv2
				call Get_param_Vr(ff, j, Vr)
				do i = 1, ff%par_nv1
					S = S + Wr * Vr * ff%DistF(i, j, k)
				end do
			end do
			S = S * 2.0 * par_pi * dWx * dWr
			if(S < 0.00001) then
				ff%Q1m(:, :, k) = 0.0_8
				ff%Q1p(k) = 0.0_8
				CYCLE
			end if
			S = 0.0

			! Для Q1p
			u = sqrt((Vx - par(2))**2 + Vr**2)
			cp = sqrt(par(3)/par(1))

			if (u / cp > 7.0) then
				uz = MK_Velosity_1(u, cp)
				ff%Q1p(k) = (uz * sig(uz))
			else
				ff%Q1p(k) = (MK_int_1(u, cp))
			end if
			
			! Для Q1m
			do i = 1, ff%par_nv1
				call Get_param_Vx(ff, i, Vx)
				do j = 1, ff%par_nv2
					call Get_param_Vr(ff, j, Vr)
					S = 0.0

					do ii = 1, ff%par_nv1
						call Get_param_Vx(ff, ii, Wx)
						do jj = 1, ff%par_nv2
							call Get_param_Vr(ff, jj, Wr)
							A = (Vx - Wx)**2 + Vr**2 + Wr**2
							B = 2.0 * Vr * Wr
							S = S + Wr * ff%DistF(ii, jj, k) * Tab1_Get(A, B)
						end do
					end do

					ff%Q1m(i, j, k) = S * 2.0_8 * dWx * dWr
				end do
			end do
		end do
		
	end subroutine Calc_Q

	subroutine Print_fx(f, x, name)
		TYPE (DistF), intent(in) :: f
		real(8), intent(in) :: x
		character(len=5), intent(in) :: name
		integer :: i, j
		real(8) :: S, Vx, Vr, ff

		open(1, file = "Dist_f_" // name // ".txt")

		do i = 1, f%par_nv1
			call Get_param_Vx(f, i, Vx)
			S = 0.0
			do j = 1, f%par_nv2
				call Get_param_Vr(f, j, Vr)
				call Get_Func(f, i, j, x, ff)
				S = S + ff * Vr
			end do
			S = S * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2
			WRITE (1, *) Vx, S, 1.0/(f%par_c * par_sqrtpi) * exp(-(Vx - f%par_Usr)**2/f%par_c**2)
		end do

		close(1)

	end subroutine Print_fx

	subroutine Print_fx_k(f, k, name, Time)
		TYPE (DistF), intent(in) :: f
		integer, intent(in) :: k
		real(8), intent(in) :: Time
		character(len=5), intent(in) :: name
		integer :: i, j
		real(8) :: S, Vx, Vr, ff, tf, x, tf2

		open(1, file = "Dist_f_k" // name // ".txt")

		write(1, *)  "TITLE = 'HP'  VARIABLES = u1, u2, u3, u4"

		do i = 1, f%par_nv1
			call Get_param_Vx(f, i, Vx)
			S = 0.0
			do j = 1, f%par_nv2
				call Get_param_Vr(f, j, Vr)
				!call Get_Func(f, i, j, x, ff)
				ff = f%DistF(i, j, k)
				S = S + ff * Vr
			end do
			S = S * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2
			tf = 1.0/(f%par_c * par_sqrtpi) * exp(-(Vx - f%par_Usr)**2/f%par_c**2)
			call Get_param_x(f, k, x)
			tf2 = 0.0
			if(Vx * Time >= x) tf2 = tf
			WRITE (1, *) Vx, S, tf, tf2
		end do

		close(1)

	end subroutine Print_fx_k

	subroutine Print_rho(f, name)
		TYPE (DistF), intent(in) :: f
		character(len=5), intent(in) :: name
		integer :: i, j, k
		real(8) :: S, Vx, Vr, ff, x

		open(1, file = "rho_" // name // ".txt")

		do k = 1, f%par_n
			S = 0.0
			call Get_param_x(f, k, x)
			do i = 1, f%par_nv1
				call Get_param_Vx(f, i, Vx)
				do j = 1, f%par_nv2
					call Get_param_Vr(f, j, Vr)
					call Get_Func(f, i, j, x, ff)
					S = S + ff * Vr
				end do
			end do
			S = S * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2 * (f%par_Rv1 - f%par_Lv1)/f%par_nv1
			WRITE (1, *) x, S
		end do

		close(1)

	end subroutine Print_rho

	subroutine Print_GD(g, name, TT)
		TYPE (GD), intent(in out) :: g
		character(len=5), intent(in) :: name
		real(8), intent(in) :: TT
		integer :: i
		real(8) :: x, c0

		c0 = sqrt(g%par_ggg * 0.3/1.0)

		open(1, file = "GD_" // name // ".txt")
		write(1, *) "TITLE = 'HP'  VARIABLES = X, RHO, U, P, Mach, Uteor"

		do i = 1, g%par_n - 1
			x = (get_coordinat_yzel(i, g%par_L, g%par_R, g) + get_coordinat_yzel(i + 1, g%par_L, g%par_R, g))/2.0
			if(x <= -2.0 * c0/(g%par_ggg - 1) * TT) then
				WRITE (1, *) x, g%par(:, i, 1),  dabs(g%par(2, i, 1))/sqrt(g%par_ggg * g%par(3, i, 1)/g%par(1, i, 1)), &
					0.0
			else if(x <= c0 * TT) then
				WRITE (1, *) x, g%par(:, i, 1),  dabs(g%par(2, i, 1))/sqrt(g%par_ggg * g%par(3, i, 1)/g%par(1, i, 1)), &
					2.0/(g%par_ggg + 1.0) * x/TT - 2.0 * c0/(g%par_ggg + 1.0)
			else
				WRITE (1, *) x, g%par(:, i, 1),  dabs(g%par(2, i, 1))/sqrt(g%par_ggg * g%par(3, i, 1)/g%par(1, i, 1)), &
					0.0
			end if
		end do

		close(1)

	end subroutine Print_GD

	real(8) pure function MK_int_1_f1(x)

		real(8), intent (in) :: x
	
		if (x <= 1.0) then
			MK_int_1_f1 =  6.283185155644284 + 0.000024846677279866114 * x + &
				2.0934329078277405 * x * x + 0.008055998193903208 * x * x * x - &
				0.2355169235647438 * x * x * x * x + 0.03820480582423355 * x * x * x * x * x + &
				0.006992274370591744 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f1 =  6.437524091973454 - 0.6331520099380095 * x + &
				3.1348881317268997 * x * x - 0.8454201478027856 * x * x * x + &
				0.1004702004260311 * x * x * x * x + 0.0009895488638964746 * x * x * x * x * x - &
				0.000920750276197054 * x * x * x * x * x * x
		else if (x <= 5) then
			MK_int_1_f1 =  4.4920780630505135 + 2.5133093267020654 * x + &
				1.1327223176567935 * x * x - 0.24648691152318875 * x * x * x + &
				0.031326738629523766 * x * x * x * x - 0.0021366031960331384 * x * x * x * x * x + &
				0.00005954097505746697 * x * x * x * x * x * x
		else if (x <= 7) then
			MK_int_1_f1 =  1.9138683588136232 + 5.350374732905213 * x - &
				0.16380205801427633 * x * x + 0.06765898334856263 * x * x * x - &
				0.011071118267864083 * x * x * x * x + 0.0008673476933852199 * x * x * x * x * x - &
				0.00002691859374483661 * x * x * x * x * x * x
		else if (x <= 50.0) then
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		else
			MK_int_1_f1 =  1.3138472469154294 + 5.336877156136497 * x + &
				0.020286308991329983 * x * x - 0.9780973533544137 * (x / 10.0)**3 + &
				0.26354051936651874 * (x / 10.0)**4 - 0.03711733070841712 * (x / 10.0)**5 + &
				0.002120935433043921 * (x / 10.0)**6
		end if
				 
		return
	end function MK_int_1_f1
	
	real(8) pure function MK_int_1_f2(x)
		real(8), intent (in) :: x
		
		if (x <= 1.0) then
			MK_int_1_f2 =  1.328216167939543 - 0.000004545681954848391 * x + &
				2.537368073155103 * x * x - 0.0020584991728545624 * x * x * x - &
				0.03742568018912792 * x * x * x * x - 0.010312136385277346 * x * x * x * x * x + &
				0.002767736179209713 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f2 =  1.2959616295629246 + 0.1533684067037866 * x + &
				2.2354849981206106 * x * x + 0.3113395567715921 * x * x * x - &
				0.21656309882941488 * x * x * x * x + 0.041957500887605075 * x * x * x * x * x - &
				0.0029978773724628604 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f2 =  1.903643456971281 - 1.4801836911099535 * x + 3.973958664572268 * x * x - &
				0.6482729779428982 * x * x * x + 0.07665007314658864 * x * x * x * x - &
				0.005369758193338703 * x * x * x * x * x + 0.00016605531871992049 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f2 =  -4.484415105552316 + 5.3747429756690766 * x + &
				0.8892806582308143 * x * x + 0.09767316152573671 * x * x * x - &
				0.025704749778475783 * x * x * x * x + 0.0021937998296249206 * x * x * x * x * x - &
				0.00006928845984076111 * x * x * x * x * x * x
		end if
				 
		return
	end function MK_int_1_f2
	
	real(8) pure function MK_int_1_f3(x)
		real(8), intent (in) :: x
 
		if (x <= 1.0) then
			MK_int_1_f3 = 1.2938345594193854 - 0.000031719847351174835 * x + &
				1.3183710041280094 * x * x - 0.014150512069488197 * x * x * x + &
				0.4226114681928129 * x * x * x * x - 0.06985750969880078 * x * x * x * x * x - &
				0.015347864048406958 * x * x * x * x * x * x
		else if (x <= 3.0) then
			MK_int_1_f3 = 0.9667460440956788 + 1.336271810704016 * x - &
				0.8687355257991665 * x * x + 1.7676868273627229 * x * x * x - &
				0.2731222764016417 * x * x * x * x + 0.004801770033831665 * x * x * x * x * x + &
				0.001780776080720323 * x * x * x * x * x * x
		else if (x <= 5.0) then
			MK_int_1_f3 = 4.760566734123174 - 5.048204299463048 * x + 3.332342585228025 * x * x + &
				0.47584339615235993 * x * x * x - 0.12072786272726124 * x * x * x * x + &
				0.011870955604980658 * x * x * x * x * x - 0.0004580199652402304 * x * x * x * x * x * x
		else if (x <= 7.0) then
			MK_int_1_f3 = 9.370493362469261 - 10.848615619383615 * x + 6.423326878282571 * x * x - &
				0.4148977656870439 * x * x * x + 0.025300923044176957 * x * x * x * x - &
				0.0010108688120876522 * x * x * x * x * x + 0.00001864423130429156 * x * x * x * x * x * x
		end if
			 
		return
	end function MK_int_1_f3
	
	real(8) pure function MK_int_1(x, cp)
		real(8), intent (in) :: x, cp
		real(8) :: b
		
		b = 1.0 - par_a_2 * log(cp)
		MK_int_1 = (cp / (par_sqrtpi**3)) * (b * b * MK_int_1_f1(x / cp) - &
			2.0 * par_a_2 * b * MK_int_1_f2(x / cp) + par_a_2**2 * MK_int_1_f3(x / cp))
	end function MK_int_1

	real(8) pure function MK_Velosity_1(u, cp)
		real(8), intent (in) :: u, cp
		
		if (u < 0.00001) then
			MK_Velosity_1 = 2.0 * cp / par_sqrtpi + 2.0 * u * u / (3.0 * cp * par_sqrtpi) - &
				u * u * u * u / (15.0 * cp * cp * cp * par_sqrtpi)
		else
			MK_Velosity_1 =  exp(-u * u / cp**2) * cp / par_sqrtpi + (u + (cp**2) / (2.0 * u)) * erf(u / cp)
		end if
		
	end function MK_Velosity_1

end module Distfunc