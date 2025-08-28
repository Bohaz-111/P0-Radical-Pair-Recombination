program SW_n16
    use iso_fortran_env
    implicit none
    
    integer, parameter :: dp = real64
    integer, parameter :: time_step = 3000
    real(dp), parameter :: t_end = 30.0_dp
    real(dp), parameter :: omega = 0.5_dp
    real(dp), parameter :: a(16) = [ -0.999985_dp, -0.7369246_dp, 0.511210_dp, -0.0826998_dp, &
            0.0655341_dp, -0.562082_dp, -0.905911_dp, 0.357729_dp, &
            0.358593_dp, 0.869386_dp, -0.232996_dp, 0.0388327_dp, &
            0.661931_dp, -0.930856_dp, -0.893077_dp, 0.0594001_dp ]
    
    real(dp) :: tau, omega_star
    real(dp), dimension(time_step) :: t, R_xx, R_xy, R_zz
    integer :: i
    
    tau = sqrt(8.0_dp / sum(a**2))
    omega_star = omega * tau
    
    t = [(i * t_end / real(time_step, dp), i = 1, time_step)]
    
    do i = 1, time_step
        R_xx(i) = R_xx_func(t(i), omega_star, tau)
        R_xy(i) = R_xy_func(t(i), omega_star, tau)
        R_zz(i) = R_zz_func(t(i), omega_star, tau)
    end do
    
    call write_data('SW_R_xx_n16.dat', t, R_xx)
    call write_data('SW_R_xy_n16.dat', t, R_xy)
    call write_data('SW_R_zz_n16.dat', t, R_zz)

contains

    function R_xx_func(t, omega_star, tau) result(R_xx)
        real(dp), intent(in) :: t, omega_star, tau
        real(dp) :: R_xx
        real(dp) :: t_star
        
        t_star = t / tau
        R_xx = (omega_star*(2.0_dp + exp(-t_star**2)*((omega_star**2-2.0_dp)*cos(omega_star*t_star) - &
                2.0_dp*omega_star*t_star*sin(omega_star*t_star))) - &
                4.0_dp*trapezoidal_rule(f1, 0.0_dp, t_star, 1000)) / (2.0_dp*omega_star**3)
    end function R_xx_func

    function R_xy_func(t, omega_star, tau) result(R_xy)
        real(dp), intent(in) :: t, omega_star, tau
        real(dp) :: R_xy
        real(dp) :: t_star
        
        t_star = t / tau
        R_xy = exp(-t_star**2) * (2.0_dp*omega_star*t_star*cos(omega_star*t_star) + &
               (omega_star**2-2.0_dp)*sin(omega_star*t_star)) / (2.0_dp*omega_star**2)
    end function R_xy_func

    function R_zz_func(t, omega_star, tau) result(R_zz)
        real(dp), intent(in) :: t, omega_star, tau
        real(dp) :: R_zz
        real(dp) :: t_star
        
        t_star = t / tau
        R_zz = (omega_star*(omega_star**2 + 4.0_dp*exp(-t_star**2)*cos(omega_star*t_star) - 4.0_dp) + &
               8.0_dp*trapezoidal_rule(f1, 0.0_dp, t_star, 1000)) / (2.0_dp*omega_star**3)
    end function R_zz_func

    function trapezoidal_rule(f, a, b, n) result(integral)
        real(dp), intent(in) :: a, b
        integer, intent(in) :: n
        real(dp) :: integral

        interface
            function f(x) result(y)
                use iso_fortran_env
                integer, parameter :: dp = real64
                real(dp), intent(in) :: x 
                real(dp) :: y 
            end function f
        end interface

        real(dp) :: h, sum, x
        integer :: i
        
        h = (b - a) / real(n, dp)
        sum = 0.5_dp * (f(a) + f(b))
        
        do i = 1, n-1
            x = a + i * h
            sum = sum + f(x)
        end do
        
        integral = h * sum
    end function trapezoidal_rule

    function f1(x) result(y)
        real(dp), intent(in) :: x 
        real(dp) :: y
        
        y = exp(-x**2) * sin(omega_star * x)
    end function f1

    subroutine write_data(filename, t, data)
        character(len=*), intent(in) :: filename
        real(dp), dimension(:), intent(in) :: t, data
        integer :: i, unit
        
        open(newunit=unit, file=filename, status='replace')
        do i = 1, size(t)
            write(unit,*) t(i), data(i)
        end do
        close(unit)
    end subroutine write_data

end program SW_n16
