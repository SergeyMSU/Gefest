
module Distfunc
	USE STORAGE
	USE Solvers
	USE My_func
	USE cgod
	USE Tab1
	USE OMP_LIB
	USE Emath
	USE MK
	implicit none
	
	contains
	
	subroutine Init_Func(f)
		! Инициализации функции распределения
		TYPE (DistF), intent(in out) :: f
		integer :: i, j, k, n
		real(8) :: Vx, Vr, x, tt
		logical:: exists
		n = 0
		allocate(f%DistF(f%par_nv1, f%par_nv2, f%par_n))



		if(ALLOCATED(QQ) == .False.) then
			allocate(QQ(f%par_nv1, f%par_nv2, f%par_n))
			allocate(Q2(f%par_n))
			allocate(Q3(f%par_n))
			QQ = 0.0
			Q2 = 0.0
			Q3 = 0.0

			allocate(Q1m(f%par_nv1, f%par_nv2, f%par_n))
			allocate(Q1p(f%par_nv1, f%par_nv2, f%par_n))

			allocate(Q1mHH(f%par_nv1, f%par_nv2, f%par_n))
			allocate(Q1pHH(f%par_nv1, f%par_nv2, f%par_n))

			Q1m = 0.0
			Q1p = 0.0

			Q1mHH = 0.0
			Q1pHH = 0.0

		end if
		
		f%DistF = 0.0_8
		do k = 1, f%par_n
			call Get_param_x(f, k, x)
			do j = 1,  f%par_nv2
				call Get_param_Vr(f, j, Vr)
				do i = 1, f%par_nv1
					call Get_param_Vx(f, i, Vx)
					if(x <= 0.0) f%DistF(i, j, k) = f%par_nH * f_maxwell(Vx, Vr, f%par_Usr, f%par_c)
				end do
			end do
		end do

		if(.False.)then!(ALLOCATED(time_step) == .False.) then
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
			g%par(1, i, 1) = pl_rho
			g%par(1, i, 2) = pl_rho

			g%par(2, i, 1) = pl_u
			g%par(2, i, 2) = pl_u

			g%par(3, i, 1) = pl_p
			g%par(3, i, 2) = pl_p
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
		real(8) :: dx, dx_do, newL, newR, w, QQ2, QQ3, x
		real(8) :: par1(size(g%par(:, 1, 1)))
		real(8) :: par2(size(g%par(:, 1, 1)))
		real(8) :: POTOK(size(g%par(:, 1, 1)))
		real(8) :: normal, dsl, dsp, dsc, vvv
		real(8) :: qqq1(9), qqq2(9), POTOK2(9), ro, p, u, ro2, p2, u2, pp, c0, Mach
		logical :: contact, lin_Mach
		integer(4) :: kdir, idgod, KOBL

		c0 = sqrt(g%par_ggg * pl_p/pl_rho)
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
		!Mach = dabs(g%par(2, 1, now))/sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now))
		vvv = (g%par(2, 1, now) - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now)))

		!if(vvv > 0.0) vvv = - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * pl_p)
		!vvv = - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * pl_p)

		! if(Mach > 7) then
		! 	vvv = (- 7.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now)))
		! end if

		newL = g%par_L  + vvv * T
		lin_Mach = .False.

		if(newL <= f1%par_L) then
			lin_Mach = .True.
			newL = f1%par_L
		end if

		! newL = g%par_L  + (g%par(2, 1, now) - 2.0/(g%par_ggg - 1.0) * sqrt(g%par_ggg * g%par(3, 1, now)/g%par(1, 1, now))) * T
		! if(newL < -1.5) newL = -1.5

		newR = g%par_R

		dx = get_coordinat_yzel(2, newL, newR, g) - get_coordinat_yzel(1, newL, newR, g)
		dx_do = get_coordinat_yzel(2, g%par_L, g%par_R, g) - get_coordinat_yzel(1, g%par_L, g%par_R, g)


		do i = 1, g%par_n - 1
			
			par1 = g%par(:, i, now)
			POTOK = 0.0


			! Найдём координаты центра ячейки
			x = (get_coordinat_yzel(i + 1, g%par_L, g%par_R, g) + get_coordinat_yzel(i, g%par_L, g%par_R, g))/2.0 
			call Get_Q(x, QQ2, QQ3)

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
				
				!if(j == 2 .and. i == 1 .and. lin_Mach == .False.) contact = .True.

				!if(.False.) then!
				if(j == 2 .and. i == 1 .and. lin_Mach == .False.) then
					POTOK2 = 0.0
				else
					POTOK2 = 0.0
					!if(j == 2 .and. i == 1 .and. lin_Mach == .True.) w = 0.0_8

					! call cgod3d(KOBL, 0, 0, 0, kdir, idgod, &
					! 			normal, 0.0_8, 0.0_8, 1.0_8, &
					! 			w, qqq1(1:8), qqq2(1:8), &
					! 			dsl, dsp, dsc, 1.0_8, 1.66666666666666_8, &
					! 			POTOK2, kontact_ = contact)
					! if (idgod == 2) then
                    !             POTOK2 = 0.0
                    !             call chlld_Q(2, normal, 0.0_8, 0.0_8, &
                    !                 w, qqq1, qqq2, dsl, dsp, dsc, POTOK2, .False.)
					! end if
					call chlld_Q(2, normal, 0.0_8, 0.0_8, &
                                    w, qqq1, qqq2, dsl, dsp, dsc, POTOK2, .False.)

					
				end if

				!$omp critical
					g%time_step = min( g%time_step,  0.2 * dx/( max(dabs(dsl), dabs(dsp)) + dabs(w) ))
				!$omp end critical

				!if (idgod == 2) print*, "ERROR 4i0u43h9h43t3r434r3"

				POTOK(1) = POTOK(1) + POTOK2(1)
				POTOK(2) = POTOK(2) + POTOK2(2)
				POTOK(3) = POTOK(3) + POTOK2(5)
			end do

			ro = par1(1)
			p = par1(3)
			u = par1(2)

			QQ2 = QQ2 * ro
			QQ3 = QQ3 * ro

			ro2 = par1(1) * dx_do/dx - T * (POTOK(1) / dx)
			if(ro2 < 0.0) then
				print*, "ro2 < 0", ro2, ro, i, lin_Mach, QQ2, QQ3
				ro2 = 0.000001
			end if
			u2 = (ro * u * dx_do/dx - T * ( POTOK(2) / dx - QQ2)) / ro2
			p2 = ((  ( p / (g%par_ggg - 1.0) + 0.5 * ro * u**2)  * dx_do/dx   &
                        - T * (POTOK(3)/ dx - QQ3) ) - 0.5 * ro2 * u2**2 ) * (g%par_ggg - 1.0)
			if(p2 < 0.0) then
				!print*, "p2 < 0", p2, p, i
				p2 = 0.000001
			end if
			
			g%par(1, i, now2) = ro2
			g%par(3, i, now2) = p2
			g%par(2, i, now2) = u2

		end do

		g%par_L = newL

		if(g%par_L < f1%par_L) print*, "ERROR LLL fwv4b643n6436"

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

	subroutine Get_Q(x, QQ1, QQ2)
		real(8), intent(out) :: QQ1, QQ2      ! Значение функции в этой точке
		real(8), intent(in) :: x

		integer :: k1, k2
		real(8) :: dd, al, Vx, Vr, L, R
		
		if(x >= f1%par_R) then
			QQ1 = 0.0_8
			QQ2 = 0.0_8
			return
		end if
		
		if(x <= f1%par_L) then
			QQ1 = 0.0_8
			QQ2 = 0.0_8
			return
		end if
		
		dd = (f1%par_R - f1%par_L)/f1%par_n
		k1 = INT((x - f1%par_L)/dd) + 1
		k2 = k1 + 1
		
		L = f1%par_L + k1 * dd - dd/2.0
		R = f1%par_L + k2 * dd - dd/2.0
		if(L > x) then
			k1 = k1 - 1
			k2 = k2 - 1
			L = f1%par_L + k1 * dd - dd/2.0
			R = f1%par_L + k2 * dd - dd/2.0
		end if
		al = (x - L)/(R - L)
		
		if(k2 > f1%par_n) then
			QQ1 = (1.0_8 - al) * Q2(k1)
			QQ2 = (1.0_8 - al) * Q3(k1)
			return
		end if
		
		if(k1 < 1) then
			QQ1 = 0.0_8
			QQ2 = 0.0_8
			return
		end if
		
		QQ1 = (1.0_8 - al) * Q2(k1) + al * Q2(k2)
		QQ2 = (1.0_8 - al) * Q3(k1) + al * Q3(k2)
		return
	end subroutine Get_Q

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

	subroutine play_GD(T, now)
		! Запускаем газовую динамику СТРОГО до времени T
		real(8), intent(in) :: T
		integer, intent(in out) :: now   ! При вхождении показывает какие параметры актуальны (соотвественно другие надо менять)
		real(8) :: all_t, dt

		all_t = 0.0

		do while (all_t < T)
			dt = gd1%time_step
			if(all_t + dt > T) dt = T - all_t
			call Start_GD(dt, gd1, now)
			all_t = all_t + dt
			now = mod(now, 2) + 1
		end do

		!print*, "dt = ", dt, T, gd1%par_L
	end subroutine play_GD
	
	subroutine Integrate_Protiv_potoka(ff1, ff2, ff3, TT, now)
		TYPE (DistF), intent(in out) :: ff1, ff2, ff3
		real(8), intent(in out) :: TT  ! Общее время решения
		integer, intent(in out) :: now  ! Общее время решения
		integer :: i, j, k, st, step, m, nnn
		integer :: nn(ff1%par_nv1)
		real(8) :: dt_global, Vx_max, dx, nu, Vx, Vr, Time, Vx_min, dt, T1, T2
		real(8) :: nu_ex, N_ex, x, dt_ex, dt_ex2
		real(8) :: dQ
		real(8) :: u2p, u2m, u1, u1p, u1m, u1mm, omeg, umm, um, u, up, upp
		logical :: metod_ytochn
 
		metod_ytochn = .False.   !! Метод Королькова аналитического уточнения функции распределения 

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
		dt_global = dx/dabs(Vx_max)                                    ! Шаг по времени (для функции распределения)
		call Get_param_Vx(ff1, 1, Vx_max)
		dt_global = min(dt_global, dx/dabs(Vx_max)) !!!!!
		dt_ex = dt_global/2.0
		!print*, "dt_ex = ", dt_ex

		call Get_param_Vx(ff1, ff1%par_nv1/2 + 1, Vx_min)
		if(metod_ytochn) dt_global = dx/Vx_min                                    ! Шаг по времени (для функции распределения)
		
		! if (metod_ytochn) then
		! 	TT = dt_global
		! else
		! 	TT = dt_global * 300
		! end if
		T1 = 0.0
		dt_ex2 = 0.0

		! Сохраняем исходную функцию распределение (чтобы от неё шагать точно)
		do i = 1, ff1%par_nv1
			do k = 3, ff1%par_n-2
				do j = 1, ff1%par_nv2
					ff3%DistF(i, j, k) = ff1%DistF(i, j, k)
				end do
			end do
		end do

		QQ = 0.0_8

		nnn = 300
		if(metod_ytochn) nnn = size(time_step)

		nnn = (INT(TT/dt_global) + 2) * 2
		dt_global = TT/nnn
		dt_ex = dt_global/2.0

		do step = 1, nnn
			if(mod(step, 1) == 0) print*, "step = ", step, "From = ", nnn

			if(metod_ytochn) then
				T2 = time_step(step) * dt_global
			else
				T2 = dt_global * step
			end if

			dt = T2 - T1
			dt_ex2 = dt_ex2 + dt
			T1 = T2
			nu_ex = 0.0_8
			!print*, "dt = ", dt

			call play_GD(dt, now)


			if(dt_ex < dt_ex2) then
				dt_ex2 = 0.0
				print*, "Start Calc_Q"
				call Calc_Q(ff1, gd1, now) !! -----------------------------------------------
				print*, "End Calc_Q"
			end if

			! Делаем сначала стандартный шаг по времени dt
			do i = 1, ff1%par_nv1
				call Get_param_Vx(ff1, i, Vx)
				nu = Vx * dt/dx

				omeg = (4.0_8 * nu**2 + 1.0_8) * (4.0_8 - nu**2)/5.0_8
				!omeg = (4 * nu**2 - nu**4)

				if(dabs(nu) > 1.0) print*, "PROBLEM 9u97oph[ijouhi"

				do k = 3, ff1%par_n-2
					do j = 1, ff1%par_nv2
						call Get_param_Vr(ff1, j, Vr)
						nu_ex = Q1p(i, j, k) + Q1pHH(i, j, k)! proton(1, k) * ff1%Q1p(i, j, k)
						N_ex = Q1m(i, j, k) + Q1mHH(i, j, k)! proton(1, k) * f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k)

						if(metod_ytochn == .False.) then
							if(.false.) then
								if(Vx > 0.0) then
									
									!dQ = dt * proton(1, k) * (-ff1%Q1p(i, j, k) * ff1%DistF(i, j, k) + &
									!						f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
									ff2%DistF(i, j, k) = ff1%DistF(i, j, k) - nu * (ff1%DistF(i, j, k) - ff1%DistF(i, j, k-1)) + &
										0.5_8 * nu * (nu - 1.0_8) * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k - 1) + ff1%DistF(i, j, k-2))! &
										!+ dQ
									!QQ(i, j, k) = QQ(i, j, k) + dQ 
									!!if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
								else
									!dQ = dt * proton(1, k) * (-ff1%Q1p(i, j, k) * ff1%DistF(i, j, k) + &
									!							f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * ff1%Q1m(i, j, k))
									ff2%DistF(i, j, k) = ff1%DistF(i, j, k) + nu * (3.0_8/2.0_8 * ff1%DistF(i, j, k) - &
										2.0_8 * ff1%DistF(i, j, k + 1) + 0.5_8 * ff1%DistF(i, j, k + 2)) + &
										0.5_8 * nu**2 * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k + 1) + ff1%DistF(i, j, k+2))! &
										!+ dQ
									!QQ(i, j, k) = QQ(i, j, k) + dQ 
									!!if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
								end if
							else
								umm = ff1%DistF(i, j, k - 2)
								um = ff1%DistF(i, j, k - 1)
								u = ff1%DistF(i, j, k)
								up = ff1%DistF(i, j, k + 1)
								upp = ff1%DistF(i, j, k + 2)

								u1mm = umm - 2.0_8/3.0_8 * nu * (um - umm)
								u1m = um - 2.0_8/3.0_8 * nu * (u - um)
								u1 = u - 2.0_8/3.0_8 * nu * (up - u)
								u1p = up - 2.0_8/3.0_8 * nu * (upp - up)

								u2m = (um + u1m - 2.0_8/3.0_8 * nu * (u1m - u1mm))/2.0_8
								u2p = (up + u1p - 2.0_8/3.0_8 * nu * (u1p - u1))/2.0_8
								ff2%DistF(i, j, k) = u - nu * (-2 * upp + 7 * up - 7 * um + 2 * umm)/24.0_8 - &
									3.0_8/8.0_8 * nu * (u2p - u2m) - omeg/24.0_8 * (upp - 4 * up + 6 * u - 4 * um + umm)
									
								if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
							end if
						else
							dQ = dt * (-nu_ex * ff1%DistF(i, j, k) + N_ex)/par_N
							! if(Vx > 0.0) then
							! 	ff2%DistF(i, j, k) = ff1%DistF(i, j, k) - nu * (ff1%DistF(i, j, k) - ff1%DistF(i, j, k-1)) + &
							! 		0.5_8 * nu * (nu - 1.0_8) * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k - 1) + ff1%DistF(i, j, k-2)) &
							! 		+ dQ
							! 	QQ(i, j, k) = QQ(i, j, k) + dQ 
							! else
							! 	ff2%DistF(i, j, k) = ff1%DistF(i, j, k) + nu * (3.0_8/2.0_8 * ff1%DistF(i, j, k) - &
							! 		2.0_8 * ff1%DistF(i, j, k + 1) + 0.5_8 * ff1%DistF(i, j, k + 2)) + &
							! 		0.5_8 * nu**2 * (ff1%DistF(i, j, k) - 2.0_8 * ff1%DistF(i, j, k + 1) + ff1%DistF(i, j, k+2)) &
							! 		+ dQ
							! 	QQ(i, j, k) = QQ(i, j, k) + dQ 
							! end if

							umm = ff1%DistF(i, j, k - 2)
							um = ff1%DistF(i, j, k - 1)
							u = ff1%DistF(i, j, k)
							up = ff1%DistF(i, j, k + 1)
							upp = ff1%DistF(i, j, k + 2)

							u1mm = umm - 2.0_8/3.0_8 * nu * (um - umm)
							u1m = um - 2.0_8/3.0_8 * nu * (u - um)
							u1 = u - 2.0_8/3.0_8 * nu * (up - u)
							u1p = up - 2.0_8/3.0_8 * nu * (upp - up)

							u2m = (um + u1m - 2.0_8/3.0_8 * nu * (u1m - u1mm))/2.0_8
							u2p = (up + u1p - 2.0_8/3.0_8 * nu * (u1p - u1))/2.0_8
							ff2%DistF(i, j, k) = u - nu * (-2 * upp + 7 * up - 7 * um + 2 * umm)/24.0_8 - &
								3.0_8/8.0_8 * nu * (u2p - u2m) - omeg/24.0_8 * (upp - 4 * up + 6 * u - 4 * um + umm) + dQ

						end if

						if(metod_ytochn == .False.) then
							if(nu_ex > 0.0000001) then
								ff2%DistF(i, j, k) = ff2%DistF(i, j, k) * exp(-dt * nu_ex) + N_ex/nu_ex * (1.0 - exp(-dt * nu_ex))
							end if
						end if

						if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
					end do ! j
				end do! k
			end do !  i

			!! Уточняем функцию для нужных i
			if(metod_ytochn) then
			!!if(.False.) then
				do st = 1, ff1%par_nv1
					i = step_algoritm(st, step)
					if(i == 0) EXIT
					call Get_param_Vx(ff1, i, Vx)
					do k = 3, ff1%par_n-2
						do j = 1, ff1%par_nv2
							if(Vx > 0.0) then
								ff2%DistF(i, j, k) = ff3%DistF(i, j, k-1)  + QQ(i, j, k)
							else
								ff2%DistF(i, j, k) = ff3%DistF(i, j, k+1) + QQ(i, j, k)
							end if
							QQ(i, j, k) = 0.0_8
							if(ff2%DistF(i, j, k) < 0.0) ff2%DistF(i, j, k) = 0.0
						end do
					end do

					do k = 3, ff1%par_n-2
						do j = 1, ff1%par_nv2
							ff3%DistF(i, j, k) = ff2%DistF(i, j, k)
						end do
					end do
				end do
			end if

			do i = 1, ff1%par_nv1
				do k = 3, ff1%par_n-2
					do j = 1, ff1%par_nv2
						ff1%DistF(i, j, k) = ff2%DistF(i, j, k)
					end do
				end do
			end do

		end do ! while

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
						+ time * proton(1, k) * (-Q1p(i, j, k) * ptr_f1%DistF(i, j, k) + &
												f_maxwell(Vx, Vr, proton(2, k), proton(3, k)) * Q1m(i, j, k))
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

	function Calc_Q1mHH(wx, wr, potok, k)
		real(8), intent(in) :: wx, wr
		integer, intent(in) :: potok, k  ! k - в какой точке по x берём значения функции распределения
		real(8) :: Calc_Q1mHH, S
		real(8) :: ksi1, ksi2, ksi3, ksi4, ksi5, eps, phi, chi, g, the
        real(8) :: wwx, wwr, wwx1, wwr1
        real(8) :: du, dv
		integer :: N, i1, j1, i, j

		S = 0.0
		du = (f1%par_Rv1 - f1%par_Lv1)/f1%par_nv1
		dv = (f1%par_Rv2 - f1%par_Lv2)/f1%par_nv2
        do N = 1, 6400
            call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi3)
            call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi5)

            eps = 2.0 * par_pi * ksi1
            phi = 2.0 * par_pi * ksi2
            the = acos(1.0 - 2.0 * ksi3)

            g = MK_g_Get(ksi4)
            chi = MK_chi_Get(ksi5, g)

            call MK_exch_velos(wx, wr, g, g * sin(the) * cos(phi), g  * sin(the) * sin(phi), g * cos(the), chi, eps, wwx, wwr, wwx1, wwr1)
            
			if(wwx < f1%par_Lv1 .or. wwx1 < f1%par_Lv1 .or. wwx > f1%par_Rv1 .or. wwx1 > f1%par_Rv1) CYCLE
			if(wwr > f1%par_Rv2 .or. wwr1 > f1%par_Rv2) CYCLE

			i = INT((wwx - f1%par_Lv1)/du)
			if(i < 1) i = 1
			if(i > f1%par_nv1) i = f1%par_nv1

			i1 = INT((wwx1 - f1%par_Lv1)/du)
			if(i1 < 1) i1 = 1
			if(i1 > f1%par_nv1) i1 = f1%par_nv1

			j = INT((wwr)/dv)
			if(j < 1) j = 1
			if(j > f1%par_nv2) j = f1%par_nv2

			j1 = INT((wwr1)/dv)
			if(j1 < 1) j1 = 1
			if(j1 > f1%par_nv2) j1 = f1%par_nv2

			S = S + f1%DistF(i, j, k) * f1%DistF(i1, j1, k)
			!S = S + f_maxwell(wwx, wwr, -0.1_8, 0.2_8) * f_maxwell(wwx1, wwr1, -0.1_8, 0.2_8) 
        end do

        S = S * (MK_norm_A/6400)
		Calc_Q1mHH = S
	end function Calc_Q1mHH

	subroutine Calc_Q1mHH_all(potok, k)
		integer, intent(in) :: potok, k  ! k - в какой точке по x берём значения функции распределения
		real(8) :: S
		real(8) :: ksi1, ksi2, ksi3, ksi4, ksi5, eps, phi, chi, g, the
        real(8) :: wwx, wwr, wwx1, wwr1, wx, wr
        real(8) :: du, dv
		integer :: N, i1, j1, i, j, ii, jj
		real(8) :: lx, ly, lz, ksiax, ksiay, ksiaz, gxy, wwy, wwz, wwy1, wwz1, gx, gy, gz

		du = (f1%par_Rv1 - f1%par_Lv1)/f1%par_nv1
		dv = (f1%par_Rv2 - f1%par_Lv2)/f1%par_nv2

        do N = 1, 100000
            call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi3)
            call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi5)

            eps = 2.0 * par_pi * ksi1
            phi = 2.0 * par_pi * ksi2
            the = acos(1.0 - 2.0 * ksi3)

            g = MK_g_Get(ksi4)
            chi = MK_chi_Get(ksi5, g)

			gx = g * sin(the) * cos(phi)
			gy = g  * sin(the) * sin(phi)
			gz = g * cos(the)
			gxy = sqrt(gx**2 + gy**2)
			lx = gx/g * cos(chi) - gx * gz/g/gxy * cos(eps) * sin(chi) + gy/gxy * sin(eps) * sin(chi)
			ly = gy/g * cos(chi) - gy * gz/g/gxy * cos(eps) * sin(chi) - gx/gxy * sin(eps) * sin(chi)
			lz = gz/g * cos(chi) + gxy/g * cos(eps) * sin(chi)
			ksiaz = 0.5 * (- gz)

			do ii = 1, f1%par_nv1
				call Get_param_Vx(f1, ii, wx)
				do jj = 1, f1%par_nv2
					call Get_param_Vr(f1, jj, wr)

					ksiax = 0.5 * (2.0 * wx - gx)
					ksiay = 0.5 * (2.0 * wr - gy)

					wwx = ksiax + 0.5 * g * lx
					wwy = ksiay + 0.5 * g * ly
					wwz = ksiaz + 0.5 * g * lz

					wwx1 = ksiax - 0.5 * g * lx
					wwy1 = ksiay - 0.5 * g * ly
					wwz1 = ksiaz - 0.5 * g * lz

					wwr = sqrt(wwy**2 + wwz**2)
					wwr1 = sqrt(wwy1**2 + wwz1**2)

					if(wwx < f1%par_Lv1 .or. wwx1 < f1%par_Lv1 .or. wwx > f1%par_Rv1 .or. wwx1 > f1%par_Rv1) CYCLE
					if(wwr > f1%par_Rv2 .or. wwr1 > f1%par_Rv2) CYCLE

					i = INT((wwx - f1%par_Lv1)/du)
					if(i < 1) i = 1
					if(i > f1%par_nv1) i = f1%par_nv1

					i1 = INT((wwx1 - f1%par_Lv1)/du)
					if(i1 < 1) i1 = 1
					if(i1 > f1%par_nv1) i1 = f1%par_nv1

					j = INT((wwr)/dv)
					if(j < 1) j = 1
					if(j > f1%par_nv2) j = f1%par_nv2

					j1 = INT((wwr1)/dv)
					if(j1 < 1) j1 = 1
					if(j1 > f1%par_nv2) j1 = f1%par_nv2

					Q1mHH(ii, jj, k) = Q1mHH(ii, jj, k) + f1%DistF(i, j, k) * f1%DistF(i1, j1, k)
				end do
			end do
        end do

        Q1mHH(:, :, k) = Q1mHH(:, :, k) * (MK_norm_A/100000)
		
	end subroutine Calc_Q1mHH_all

	subroutine Calc_Q1mHH_all_k(nH_mas)
		!integer, intent(in) :: potok  ! k - в какой точке по x берём значения функции распределения
		real(8), intent(in) :: nH_mas(f1%par_n)
		real(8) :: ksi1, ksi2, ksi3, ksi4, ksi5, eps, phi, chi, g, the
        real(8) :: wwx, wwr, wwx1, wwr1, wx, wr
        real(8) :: du, dv
		integer :: N, i1, j1, i, j, ii, jj, k, potok
		integer  :: iijj(120, 35, 120, 35) = 0
		real(8) :: lx, ly, lz, ksiax, ksiay, ksiaz, gxy, wwy, wwz, wwy1, wwz1, gx, gy, gz

		iijj = 0

		du = (f1%par_Rv1 - f1%par_Lv1)/f1%par_nv1
		dv = (f1%par_Rv2 - f1%par_Lv2)/f1%par_nv2
		potok = 1

		! !$omp parallel

		! !$omp do private(lx, ly, lz, ksiax, ksiay, ksiaz, gxy, wwy, wwz, wwy1, wwz1, gx, gy, gz, i1, j1, i, j, ii, jj, k, potok, ksi1, ksi2, ksi3, ksi4, ksi5, eps, phi, chi, g, the, wwx, wwr, wwx1, wwr1, wx, wr)
        do N = 1, 10000
			potok = (omp_get_thread_num() + 1)
            call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi1)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi2)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi3)
            call M_K_rand(sensor(1, 1, potok), sensor(2, 1, potok), sensor(3, 1, potok), ksi4)
            call M_K_rand(sensor(1, 2, potok), sensor(2, 2, potok), sensor(3, 2, potok), ksi5)

            eps = 2.0 * par_pi * ksi1
            phi = 2.0 * par_pi * ksi2
            the = acos(1.0 - 2.0 * ksi3)

            g = MK_g_Get(ksi4)
            chi = MK_chi_Get(ksi5, g)

			gx = g * sin(the) * cos(phi)
			gy = g  * sin(the) * sin(phi)
			gz = g * cos(the)
			gxy = sqrt(gx**2 + gy**2)
			lx = gx/g * cos(chi) - gx * gz/g/gxy * cos(eps) * sin(chi) + gy/gxy * sin(eps) * sin(chi)
			ly = gy/g * cos(chi) - gy * gz/g/gxy * cos(eps) * sin(chi) - gx/gxy * sin(eps) * sin(chi)
			lz = gz/g * cos(chi) + gxy/g * cos(eps) * sin(chi)
			ksiaz = 0.5 * (- gz)

			do ii = 1, f1%par_nv1
				call Get_param_Vx(f1, ii, wx)
				do jj = 1, f1%par_nv2
					call Get_param_Vr(f1, jj, wr)

					ksiax = 0.5 * (2.0 * wx - gx)
					ksiay = 0.5 * (2.0 * wr - gy)

					wwx = ksiax + 0.5 * g * lx
					wwy = ksiay + 0.5 * g * ly
					wwz = ksiaz + 0.5 * g * lz

					wwx1 = ksiax - 0.5 * g * lx
					wwy1 = ksiay - 0.5 * g * ly
					wwz1 = ksiaz - 0.5 * g * lz

					wwr = sqrt(wwy**2 + wwz**2)
					wwr1 = sqrt(wwy1**2 + wwz1**2)

					if(wwx < f1%par_Lv1 .or. wwx1 < f1%par_Lv1 .or. wwx > f1%par_Rv1 .or. wwx1 > f1%par_Rv1) CYCLE
					if(wwr > f1%par_Rv2 .or. wwr1 > f1%par_Rv2) CYCLE

					i = INT((wwx - f1%par_Lv1)/du)
					if(i < 1) i = 1
					if(i > f1%par_nv1) i = f1%par_nv1

					i1 = INT((wwx1 - f1%par_Lv1)/du)
					if(i1 < 1) i1 = 1
					if(i1 > f1%par_nv1) i1 = f1%par_nv1

					j = INT((wwr)/dv)
					if(j < 1) j = 1
					if(j > f1%par_nv2) j = f1%par_nv2

					j1 = INT((wwr1)/dv)
					if(j1 < 1) j1 = 1
					if(j1 > f1%par_nv2) j1 = f1%par_nv2

					do k = 1, f1%par_n
						if(nH_mas(k) < 0.00001) CYCLE

						! !$omp critical
						!Q1mHH(ii, jj, k) = Q1mHH(ii, jj, k) + f1%DistF(i, j, k) * f1%DistF(i1, j1, k)
						iijj(i, j, i1, j1) = iijj(i, j, i1, j1) + 1
						! !$omp end critical

					end do
				end do
			end do
        end do
		! !$omp end do
		! !$omp end parallel


        Q1mHH(:, :, :) = Q1mHH(:, :, :) * (MK_norm_A/10000)
		
	end subroutine Calc_Q1mHH_all_k

	subroutine Calc_Q(ff, g, now)
		TYPE (DistF), intent(in out), TARGET :: ff
		TYPE (GD), intent(in) :: g
		integer(4), intent(in) :: now
		integer(4) :: k, i, j, ii, jj, potok
		real(8) :: x, Vx, Vr, Wx, Wr
		real(8) :: dWx, dWr, A, B, r, ecint, nH
		real(8) :: par(3)
		real(8) :: u, S, cp, uz, SS, tab, u_proton, rho, SQ2, SQ3, SS2, tab2
		real(8) :: proton(3, ff%par_n)
		real(8) :: nH_mas(ff%par_n)
		logical :: nullf(ff%par_nv1)   !! Попробуем ввести небольшое ускорение счёта интеграллов
		real(8) :: start_time, end_time

		dWx = (ff%par_Rv1 - ff%par_Lv1)/(ff%par_nv1)
		dWr = (ff%par_Rv2 - ff%par_Lv2)/(ff%par_nv2)
		start_time = omp_get_wtime()

		Q1m = 0.0_8
		Q1p = 0.0_8
		nH_mas = 0.0_8

		Q1mHH = 0.0_8
		Q1pHH = 0.0_8

		Q2 = 0.0_8
		Q3 = 0.0_8

		call omp_set_num_threads(MK_n_potok)

		


		!$omp parallel

		!$omp do private(S, i, j, Vr, tab)
		do k = 1, ff%par_n	 !! Считаем концентрацию функции распределения в каждой точке
			S = 0.0
			do j = 1, ff%par_nv2
				call Get_param_Vr(ff, j, Vr)
				do i = 1, ff%par_nv1
					tab = ff%DistF(i, j, k)
					S = S + Vr * tab
				end do
			end do
			S = S * 2.0 * par_pi * dWx * dWr
			nH_mas(k) = S
		end do
		!$omp end do

		!$omp do private(x, par) schedule(dynamic, 2)
		do k = 1, ff%par_n	   !! Считываем параметры плазмы в каждой точке
			call Get_param_x(ff, k, x)
			call Get_GD(x, gd1, par, now)
			proton(1, k) = par(1)
			proton(2, k) = par(2)
			if(par(1) > 0.000001) then
				proton(3, k) = sqrt(par(3)/par(1))
			else
				proton(3, k) = 1.0_8
			end if
		end do
		!$omp end do
		
        !! Учёт перезарядки
		if(.True.) then
		!$omp do private(tab2, SS2, potok, nullf, i, j, ii, jj, x, Vx, Vr, Wx, Wr, A, B, u, S, cp, uz, r, ecint, SS, tab, nH, u_proton, rho, SQ2, SQ3) schedule(dynamic, 2)
		do k = 1, ff%par_n  ! Пробегаемся по пространству
		
			potok = (omp_get_thread_num() + 1)
			cp = proton(3, k)
			u_proton = proton(2, k)
			rho = proton(1, k)

			if(rho < 0.00001) then
				Q1m(:, :, k) = 0.0_8
				Q1p(:, :, k) = 0.0_8
				Q2(k) = 0.0_8
				Q3(k) = 0.0_8
				CYCLE
			end if

			nH = nH_mas(k)
			if(nH < 0.00001) then
				Q1m(:, :, k) = 0.0_8
				Q1p(:, :, k) = 0.0_8
				Q2(k) = 0.0_8
				Q3(k) = 0.0_8
				CYCLE
			end if
			
			S = 0.0

			! Для Q1p
			

			
			! Для Q1m
			SQ2 = 0.0
			SQ3 = 0.0

			do i = 1, ff%par_nv1
				!!if(nullf(i) == .False.) CYCLE
				call Get_param_Vx(ff, i, Vx)
				do j = 1, ff%par_nv2
					call Get_param_Vr(ff, j, Vr)
					S = 0.0
					SS = 0.0
					SS2 = 0.0

					do ii = 1, ff%par_nv1
						call Get_param_Vx(ff, ii, Wx)
						do jj = 1, ff%par_nv2
							call Get_param_Vr(ff, jj, Wr)
							tab = Omega1(i, j, ii, jj)! Tab1_Get(A, B)
							!!if(nullf(ii) == .True.) 
							S = S + Wr * ff%DistF(ii, jj, k) * tab
							SS = SS + Wr * f_maxwell(Wx, Wr, u_proton, cp) * tab
							!!if(nullf(i) == .True.) then
								SQ2 = SQ2 + tab * Wr * Vr * (Vx - Wx) * f_maxwell(Wx, Wr, u_proton, cp) * ff%DistF(i, j, k)
								SQ3 = SQ3 + tab * Wr * Vr * (Vx**2 + Vr**2 - Wx**2 - Wr**2) * f_maxwell(Wx, Wr, u_proton, cp) * ff%DistF(i, j, k)
							!!end if
						end do
					end do

					Q1m(i, j, k) = S * 2.0_8 * dWx * dWr! * 2.0 * par_pi
					Q1p(i, j, k) = SS * 2.0_8 * dWx * dWr!
					!! Предлагается здесь сразу умножить источники на то что надо
					Q1m(i, j, k) = Q1m(i, j, k) * rho * f_maxwell(Vx, Vr, u_proton, cp)
				end do
			end do

			Q2(k) = QKnHp * nH * SQ2 * 4.0_8 * par_pi * dWx * dWr * dWx * dWr  ! На rho надо умножить в программе газовой динамики
			Q3(k) = QKnHp * nH * SQ3 * 2.0_8 * par_pi * dWx * dWr * dWx * dWr  ! На rho надо умножить в программе газовой динамики

			!! Предлагается здесь сразу умножить источники на то что надо (чтобы в основной программе не умножать и не вычислять параметры плазмы)
			Q1p(:, :, k) = Q1p(:, :, k) * rho * KnHp
			Q1m(:, :, k) = Q1m(:, :, k) * KnHp
		end do
		!$omp end do
		end if

		!$omp end parallel

		!! Учёт HH-столкновений ----------------------------------------------------------

		if(.True.) then
			! Для Q1pHH
			if(.True.) then
			do i = 1, ff%par_nv1
				call Get_param_Vx(ff, i, Vx)
				do j = 1, ff%par_nv2
					call Get_param_Vr(ff, j, Vr)
					Q1pHH(i, j, :) = 0.0
					do ii = 1, ff%par_nv1
						call Get_param_Vx(ff, ii, Wx)
						do jj = 1, ff%par_nv2
							call Get_param_Vr(ff, jj, Wr)
							tab2 = Omega2(i, j, ii, jj)


							do k = 1, ff%par_n  ! Пробегаемся по пространству
								if(nH_mas(k) < 0.00001) then
									CYCLE
								end if
								Q1pHH(i, j, k) = Q1pHH(i, j, k) + Wr * ff%DistF(ii, jj, k) * tab2 * 2.0_8 * dWx * dWr
							end do
						end do
					end do
					!Q1pHH(i, j, k) = SS2 * 2.0_8 * dWx * dWr!
				end do
			end do
			end if

			! Для Q1mHH
			! if(.False.) then
			! do i = 1, ff%par_nv1
			! 	call Get_param_Vx(ff, i, Vx)
			! 	do j = 1, ff%par_nv2
			! 		call Get_param_Vr(ff, j, Vr)
			! 		Q1mHH(i, j, k) = Calc_Q1mHH(Vx, Vr, potok, k)
			! 	end do
			! end do
			! end if

			! Q1mHH(:, :, k) = 0.0
			! call Calc_Q1mHH_all(potok, k)
			! if(k == ff%par_n/2) then
			! 	print*, "Integr = ", Q1mHH(60, 3, k)
			! end if

			!! Предлагается здесь сразу умножить источники на то что надо (чтобы в основной программе не умножать и не вычислять параметры плазмы)
			!Q1pHH(:, :, k) = Q1pHH(:, :, k) * nH * KnHH
			!Q1mHH(:, :, k) = Q1mHH(:, :, k) * nH * KnHH
		
		
		end if

		
		call Calc_Q1mHH_all_k(nH_mas)
		print*, "Integr = ", Q1mHH(60, 3, ff%par_n/2)
		do k = 1, ff%par_n
			Q1mHH(:, :, k) = Q1mHH(:, :, k) * nH_mas(k) * KnHH
		end do

		do k = 1, ff%par_n  ! Пробегаемся по пространству
			Q1pHH(:, :, k) = Q1pHH(:, :, k) * nH_mas(k) * KnHH
			Q1mHH(:, :, k) = Q1mHH(:, :, k) * nH_mas(k) * KnHH
		end do



		!!
		! Q1pHH = 0.0
		! Q1mHH = 0.0

		end_time = omp_get_wtime()
		print *, "Time work: ", (end_time-start_time), "   in secunds"

	end subroutine Calc_Q

	subroutine Print_fx(f, x, name, Time_, ii_)
		TYPE (DistF), intent(in) :: f
		real(8), intent(in) :: x
		real(8), intent(in), OPTIONAL :: Time_
		integer, intent(in), OPTIONAL :: ii_
		character(len=5), intent(in) :: name
		character(len=5)  :: name2
		integer :: i, j, nn
		real(8) :: S, Vx, Vr, ff, cp, Time, fa

		Time = 0.0
		nn = 0

		if(PRESENT(Time_)) Time = Time_
		if(PRESENT(Time_)) nn = ii_

		write(unit=name2,fmt='(i5.5)') nn

		open(1, file = name2 //"_Dist_f_" // name // ".txt")
		write(1, *)  "TITLE = 'HP'  VARIABLES = Vx, f, fH, fp, f_analitic"

		cp = sqrt(pl_p/pl_rho)

		do i = 1, f%par_nv1
			call Get_param_Vx(f, i, Vx)
			S = 0.0
			do j = 1, f%par_nv2
				call Get_param_Vr(f, j, Vr)
				call Get_Func(f, i, j, x, ff)
				S = S + ff * Vr
			end do
			S = S * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2
			fa = 0.0
			if(x - Vx * Time < 0.0) fa = 1.0/(f%par_c * par_sqrtpi) * exp(-(Vx - f%par_Usr)**2/f%par_c**2)
			WRITE (1, *) Vx, S, f%par_nH/(f%par_c * par_sqrtpi) * exp(-(Vx - f%par_Usr)**2/f%par_c**2), &
				pl_rho/(cp * par_sqrtpi) * exp(-(Vx - pl_u)**2/cp**2), f%par_nH * fa
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

		write(1, *)  "TITLE = 'HP'  VARIABLES = u1, u2, fH, u4, fp"

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
			WRITE (1, *) Vx, S, tf, tf2, &
				1.0/(sqrt(0.3) * par_sqrtpi) * exp(-(Vx)**2/0.3)
		end do

		close(1)

	end subroutine Print_fx_k

	subroutine Print_rho(f, name)
		TYPE (DistF), intent(in) :: f
		character(len=5), intent(in) :: name
		integer :: i, j, k
		real(8) :: S, Vx, Vr, ff, x, SS, SSS, Temp

		open(1, file = "rho_" // name // ".txt")
		write(1, *)  "TITLE = 'HP'  VARIABLES = x, nH, Vx, TH"

		do k = 1, f%par_n
			S = 0.0
			SS = 0.0
			SSS = 0.0
			call Get_param_x(f, k, x)
			do i = 1, f%par_nv1
				call Get_param_Vx(f, i, Vx)
				do j = 1, f%par_nv2
					call Get_param_Vr(f, j, Vr)
					call Get_Func(f, i, j, x, ff)
					S = S + ff * Vr
					SS = SS + ff * Vr * Vx
				end do
			end do

			S = S * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2 * (f%par_Rv1 - f%par_Lv1)/f%par_nv1
			SS = SS * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2 * (f%par_Rv1 - f%par_Lv1)/f%par_nv1/S

			do i = 1, f%par_nv1
				call Get_param_Vx(f, i, Vx)
				do j = 1, f%par_nv2
					call Get_param_Vr(f, j, Vr)
					call Get_Func(f, i, j, x, ff)
					SSS = SSS + ff * Vr * ( (Vx - SS)**2 + Vr**2)
				end do
			end do
			SSS = SSS * 2.0 * par_pi * (f%par_Rv2 - f%par_Lv2)/f%par_nv2 * (f%par_Rv1 - f%par_Lv1)/f%par_nv1/3.0/S

			if(S > 0.00001) then
				WRITE (1, *) x, S, SS, SSS
			else
				WRITE (1, *) x, S, 0.0, 0.0
			end if
		end do

		close(1)

	end subroutine Print_rho

	subroutine Print_GD(g, name, TT)
		TYPE (GD), intent(in out) :: g
		character(len=5), intent(in) :: name
		real(8), intent(in) :: TT
		integer :: i
		real(8) :: x, c0, Temp

		c0 = sqrt(g%par_ggg * pl_p/pl_rho)

		open(1, file = "GD_" // name // ".txt")
		write(1, *) "TITLE = 'HP'  VARIABLES = X, RHO, U, P, Mach, Uteor, T"

		do i = 1, g%par_n - 1
			x = (get_coordinat_yzel(i, g%par_L, g%par_R, g) + get_coordinat_yzel(i + 1, g%par_L, g%par_R, g))/2.0
			Temp = 0.0
			if(x <= -2.0 * c0/(g%par_ggg - 1) * TT) then
				WRITE (1, *) x, g%par(:, i, 1),  dabs(g%par(2, i, 1))/sqrt(g%par_ggg * g%par(3, i, 1)/g%par(1, i, 1)), &
					0.0, g%par(3, i, 1)/g%par(1, i, 1)
			else if(x <= c0 * TT) then
				WRITE (1, *) x, g%par(:, i, 1),  dabs(g%par(2, i, 1))/sqrt(g%par_ggg * g%par(3, i, 1)/g%par(1, i, 1)), &
					2.0/(g%par_ggg + 1.0) * x/TT - 2.0 * c0/(g%par_ggg + 1.0), g%par(3, i, 1)/g%par(1, i, 1)
			else
				WRITE (1, *) x, g%par(:, i, 1),  dabs(g%par(2, i, 1))/sqrt(g%par_ggg * g%par(3, i, 1)/g%par(1, i, 1)), &
					0.0, g%par(3, i, 1)/g%par(1, i, 1)
			end if
		end do

		close(1)

	end subroutine Print_GD

	subroutine Save_setka_bin(num)
		! Variables
		integer, intent(in) :: num
		character(len=5) :: name
		integer :: i

		
		write(unit=name,fmt='(i5.5)') num
		
		open(1, file = "save_" // name // ".bin", FORM = 'BINARY')

		call Print_rho(f1, name)
		call Print_GD(gd1, name, 1.0_8)
		
		write(1)  f1%par_n, f1%par_nv1, f1%par_nv2, f1%par_L, f1%par_R
		write(1)  f1%DistF

		write(1) gd1%par_n
		write(1) gd1%par

		close(1)

	end subroutine Save_setka_bin

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