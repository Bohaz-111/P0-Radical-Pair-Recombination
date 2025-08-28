program SC_n16
    use iso_fortran_env
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: n_spins = 16, n_steps = 3000, n_samples = 40000
    real(dp), parameter :: pi = 3.1415926535897932384626_dp, B = 0.5_dp, t_end = 30.0_dp, dt = t_end/n_steps

    real(dp), dimension(n_spins) :: a
    real(dp), dimension(n_samples, n_steps) :: Rxx, Rxy, Rzz
    integer :: s

    call random_seed()
    a = [-0.999985_dp, -0.7369246_dp, 0.511210_dp, -0.0826998_dp, 0.0655341_dp, -0.562082_dp, -0.905911_dp, 0.357729_dp, &
    0.358593_dp, 0.869386_dp, -0.232996_dp, 0.0388327_dp, 0.661931_dp, -0.930856_dp, -0.893077_dp, 0.0594001_dp]

    do s = 1, n_samples
        call run_single_trajectory(a, B, dt, n_steps, Rxx(s,:), Rxy(s,:), Rzz(s,:))
    end do

    call save_results(n_steps, dt, Rxx, Rxy, Rzz, n_samples)

contains

    subroutine run_single_trajectory(a, B, dt, n_steps, Rxx, Rxy, Rzz)
        real(dp), intent(in) :: a(n_spins), B, dt
        integer, intent(in) :: n_steps
        real(dp), intent(out) :: Rxx(n_steps), Rxy(n_steps), Rzz(n_steps)

        real(dp), dimension(3) :: S, S0
        real(dp), dimension(n_spins, 3) :: I, I0
        integer :: idx, t_idx

        call random_on_sphere(S0)
        do idx = 1, n_spins
            call random_on_sphere(I0(idx,:))
        end do
        S0 = (sqrt(3.0_dp)/2.0_dp)*S0
        I0 = (sqrt(3.0_dp)/2.0_dp)*I0

        S = S0
        I = I0

        do t_idx = 1, n_steps
            call evolve(S, I, a, B, dt)
            Rxx(t_idx) = S0(1) * S(1)
            Rxy(t_idx) = S0(1) * S(2)
            Rzz(t_idx) = S0(3) * S(3)
        end do
    end subroutine

    subroutine evolve(S, I, a, B, dt)
        real(dp), intent(inout) :: S(3)
        real(dp), intent(inout) :: I(n_spins, 3)
        real(dp), intent(in) :: a(n_spins), B, dt
        real(dp), dimension(3) :: omega
        integer :: k,l,m

        omega = [0.0_dp, 0.0_dp, B]
        do k = 1, n_spins
            omega = omega + a(k) * I(k,:)
        end do                              ! calculate resultant nuclear spin - field vector
        call rotate(S, omega, 0.5_dp * dt) ! rotate S about omega for δt/2

        do l = 1, n_spins
            omega = a(l) * S
            call rotate(I(l,:), omega, dt) ! rotate individual I about S for δt
        end do

        omega = [0.0_dp, 0.0_dp, B]
        do m = 1, n_spins
            omega = omega + a(m) * I(m,:)
        end do
        call rotate(S, omega, 0.5_dp * dt) ! rotate S about omega for δt/2
    end subroutine

    subroutine rotate(v, axis, t)
        implicit none
        real(dp), dimension(3), intent(inout) :: v 
        real(dp), dimension(3), intent(in) :: axis
        real(dp), intent(in) :: t
        real(dp) :: ax_norm
        real(dp), dimension(3) :: unit_ax, v_para, v_perp, v_cross, v_rot

        ax_norm = sqrt(dot_product(axis, axis)) 
        unit_ax = axis / ax_norm

        v_para = dot_product(unit_ax, v) * unit_ax
        v_perp = v - v_para
        v_cross = cross_product(unit_ax, v)

        v_rot = v_para + v_perp*cos(ax_norm * t) + v_cross*sin(ax_norm * t)
        v = v_rot
    end subroutine

    function cross_product(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        
        c(1) = a(2) * b(3) - a(3) * b(2)
        c(2) = a(3) * b(1) - a(1) * b(3)
        c(3) = a(1) * b(2) - a(2) * b(1)
    end function cross_product

    subroutine random_on_sphere(v)
        real(dp), intent(out) :: v(3)
        real(dp) :: theta, phi, costheta

        call random_number(phi)
        call random_number(costheta)
        costheta = 2*costheta - 1
        theta = acos(costheta)
        phi = 2.0_dp * pi * phi

        v(1) = sin(theta) * cos(phi)
        v(2) = sin(theta) * sin(phi)
        v(3) = cos(theta)
    end subroutine

    subroutine save_results(n_steps, dt, Rxx, Rxy, Rzz, n_samples)
        integer, intent(in) :: n_steps, n_samples
        real(dp), intent(in) :: dt
        real(dp), intent(in) :: Rxx(n_samples, n_steps), Rxy(n_samples, n_steps), Rzz(n_samples, n_steps)
        integer :: t, s
        real(dp) :: time, avg_xx, avg_xy, avg_zz

        open(unit=11, file='SC_Rxx_n16.dat', status='replace')
        open(unit=12, file='SC_Rxy_n16.dat', status='replace')
        open(unit=13, file='SC_Rzz_n16.dat', status='replace')

        do t = 1, n_steps
            avg_xx = 0.0_dp
            avg_xy = 0.0_dp
            avg_zz = 0.0_dp
            do s = 1, n_samples
                avg_xx = avg_xx + Rxx(s,t)
                avg_xy = avg_xy + Rxy(s,t)
                avg_zz = avg_zz + Rzz(s,t)
            end do
            avg_xx = 2.0_dp*avg_xx / n_samples
            avg_xy = 2.0_dp*avg_xy / n_samples
            avg_zz = 2.0_dp*avg_zz / n_samples
            time = (t-1)*dt
            write(11,'(4f15.6)') time, avg_xx
            write(12,'(4f15.6)') time, avg_xy
            write(13,'(4f15.6)') time, avg_zz

        end do
        close(11)
        close(12)
        close(13)
    end subroutine

end program SC_n16
