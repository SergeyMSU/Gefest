
	
module STORAGE
	implicit none 
	
	real(8), parameter :: par_pi = acos(-1.0_8) 
    real(8), parameter :: par_sqrtpi = sqrt(par_pi)
	real(8), parameter :: par_a_2 = 0.0991561_8 
	real(8), parameter :: par_N = 20_8 
	real(8), parameter :: RRR = 5.0_8
	real(8), parameter :: par_g_max = sqrt(2.0**2 + 4.0 * 2.0**2)    ! Ограничение относительной скорости
	real(8), parameter :: KnHH = 0.0125676    ! Ограничение относительной скорости
	real(8), parameter :: KnHp = 5.14822    ! Ограничение относительной скорости
	real(8), parameter :: QKnHp = 5.14822    ! Ограничение относительной скорости

	integer, PARAMETER :: MK_n_potok = 24       ! Количество потоков для массива с датчиками случайных чисел


	integer(4) :: par_ch1(36) = 0

	real(8), allocatable :: time_step(:)
	integer(4), allocatable :: step_algoritm(:, :)  ! (f%par_nv1, n)

	real(8), allocatable :: QQ(:, :, :)   ! приращение функции распределения между шагами

	real(8), allocatable :: Q2(:)         ! Источники по перезарядке для плазмы
	real(8), allocatable :: Q3(:)

	real(8), allocatable :: Q1m(:, :, :)   ! Где сложный интегралл с функцией распределения водорода
	real(8), allocatable :: Q1p(:, :, :)   ! Где в интеграле фунция распределения протонов

	real(8), allocatable :: Q1mHH(:, :, :)   ! Где сложный интегралл с функцией распределения водорода
	real(8), allocatable :: Q1pHH(:, :, :)   ! Где в интеграле фунция распределения протонов


	real(8) :: pl_rho = 1.0
	real(8) :: pl_u = 0.0
	real(8) :: pl_p = 1.0

	
	TYPE DistF 
		integer(4) :: par_n = 700

		integer(4) :: par_nv1 = 120   !! ДОЛЖНО БЫТЬ ЧЁТНЫМ  (если поменяем, нужно менять файл time_step)
		integer(4) :: par_nv2 = 35
	
		real(8) :: par_L = -5.0_8
		real(8) :: par_R = RRR
		real(8) :: par_Lv1 = -3.0_8! -3.2_8
		real(8) :: par_Rv1 = 3.0_8
		real(8) :: par_Lv2 = 0.0_8
		real(8) :: par_Rv2 = 3.0_8
		
		real(8) :: par_Usr = 0.0_8! 2.54351_8
		real(8) :: par_c = 1.41421_8
		real(8) :: par_nH = 1.0
	
		real(8), allocatable :: DistF(:, :, :)   ! Частицы (par_nv1, par_nv2, par_n) = V1, V2, X
	END TYPE DistF


	TYPE GD 
		integer(4) :: par_n = 3000
	
		real(8) :: par_L = 0.0_8
		real(8) :: par_R = 12.0_8
		real(8) :: par_ggg = (5.0_8/3.0_8)
		real(8) :: time_step
		logical :: start = .True.
		
	
		! real(8), allocatable :: X(:)   ! Координаты узлов (границ ячеек)
		real(8), allocatable :: par(:, :, :)   ! (3, : par_n - 1, 2 t)
		!(rho, u, p)
	END TYPE GD
	
	TYPE (DistF):: f1
    TYPE (DistF):: f2
    TYPE (DistF):: f3
    TYPE (GD):: gd1

	contains 

	function sig(x)
		implicit none
		real(8), intent(in) :: x
		real(8) :: sig
		
		sig = (1.0 - par_a_2 * log(x))**2
		!sig = (1.0 - 0.135838_8 * log(x/1.33301_8))**2
		!sig = 1.0_8
	end function sig

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

	function f_maxwell(x, y, u, cp)
		implicit none
		real(8), intent(in) :: x, y, u, cp
		real(8) :: f_maxwell
		
		f_maxwell = 1.0/(cp**3 * par_pi**1.5) * exp(-(x - u)**2/cp**2 - (y)**2/cp**2)
	end function f_maxwell

	
	
end module STORAGE
	
	
	
	