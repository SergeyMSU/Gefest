
	
module STORAGE
	implicit none 
	
	real(8), parameter :: par_pi = acos(-1.0_8) 
    real(8), parameter :: par_sqrtpi = sqrt(par_pi)
	real(8), parameter :: par_a_2 = 0.13043_8 

	integer(4) :: par_ch1(36) = 0

	real(8), allocatable :: time_step(:)
	integer(4), allocatable :: step_algoritm(:, :)
	
	
	TYPE DistF 
		integer(4) :: par_n = 600
		integer(4) :: par_nv1 = 36   !! ДОЛЖНО БЫТЬ ЧЁТНЫМ
		integer(4) :: par_nv2 = 36
	
		real(8) :: par_L = -1.0_8
		real(8) :: par_R = 1.0_8
		real(8) :: par_Lv1 = -2.3_8
		real(8) :: par_Rv1 = 2.3_8
		real(8) :: par_Lv2 = 0.0_8
		real(8) :: par_Rv2 = 2.3_8
		
		real(8) :: par_Usr = 0.0_8
		real(8) :: par_c = 1.0_8
	
		real(8), allocatable :: DistF(:, :, :)   ! Частицы (par_nv1, par_nv2, par_n) = V1, V2, X
		real(8), allocatable :: Q1m(:, :, :)   ! Где сложный интегралл с функцией распределения водорода
		real(8), allocatable :: Q1p(:)   ! Где в интеграле фунция распределения протонов
	END TYPE DistF


	TYPE GD 
		integer(4) :: par_n = 2000
	
		real(8) :: par_L = 0.0_8
		real(8) :: par_R = 1.0_8
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
	end function sig

	
	
end module STORAGE
	
	
	
	