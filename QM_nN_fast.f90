program QM_nN_fast
  use iso_fortran_env
  use omp_lib
  implicit none

  integer, parameter :: dp = real64

  integer,  parameter :: n_spins = 16
  integer,  parameter :: D = 2**(n_spins + 1)       
  integer,  parameter :: NSTEPS = 300
  real(dp), parameter :: T_END  = 30.0_dp
  real(dp), parameter :: dt     = T_END / NSTEPS
  real(dp), parameter :: omega  = 0.0_dp
  integer,  parameter :: m_stages = 4              
  real(dp), parameter :: acoef(m_stages) = [ &
       0.5153528374311229_dp, -0.08578201941297365_dp, &
       0.4415830236164665_dp,  0.1288461583653842_dp ]
  real(dp), parameter :: bcoef(m_stages) = [ &
       0.1344961992774311_dp, -0.2248198030794208_dp, &
       0.7563200005156683_dp,  0.3340036032863214_dp ]

  integer,  parameter :: Z = 2**n_spins           
  real(dp) :: t(NSTEPS), Rxx(NSTEPS), Rxy(NSTEPS), Rzz(NSTEPS)
  real(dp) :: a(n_spins)
  integer  :: i, M, M_prime
  real(dp) :: t0, t1

  do i = 1, NSTEPS
    t(i) = (i-1)*dt
  end do
  Rxx = 0.0_dp;  Rxy = 0.0_dp;  Rzz = 0.0_dp

  call init_hyperfine(a) 
  call omp_set_num_threads(5)

  t0 = omp_get_wtime()

  !$OMP PARALLEL DEFAULT(none) &
  !$OMP & SHARED(a,t) &
  !$OMP & PRIVATE(M,i,M_prime) &
  !$OMP & REDUCTION(+:Rxx,Rxy,Rzz)
  block
    real(dp), allocatable :: q_u(:), p_u(:), q_d(:), p_d(:)
    real(dp), allocatable :: qtmp(:), ptmp(:), Hq(:), Hp(:)
    complex(dp) :: s1, s2, uu, ud, du, dd
    real(dp)    :: n_uu, n_ud, n_du, n_dd
    integer :: j

    allocate(q_u(D), p_u(D), q_d(D), p_d(D), qtmp(D), ptmp(D), Hq(D), Hp(D))

    !$OMP DO SCHEDULE(DYNAMIC, 8)
    do M = 1, Z
      !$OMP CRITICAL
      if (mod(M, max(1, Z/1000)) == 0 .or. M == Z) then
        print '(A,I0,A,I0,A)', 'Processing nuclear state ', M, '/', Z, '...'
      end if
      !$OMP END CRITICAL
      
      call init_ket(M, .true.,  q_u, p_u)  
      call init_ket(M, .false., q_d, p_d)  

      do i = 1, NSTEPS
        call propag_4stage(q_u, p_u, dt, acoef, bcoef, a, omega, qtmp, ptmp, Hq, Hp)
        call propag_4stage(q_d, p_d, dt, acoef, bcoef, a, omega, qtmp, ptmp, Hq, Hp)

        s1 = (0.0_dp, 0.0_dp); s2 = (0.0_dp, 0.0_dp)
        n_uu = 0.0_dp; n_ud = 0.0_dp; n_du = 0.0_dp; n_dd = 0.0_dp

        do M_prime = 1, Z
          if (abs(popcnt(M_prime-1) - popcnt(M-1)) > 1) cycle

          uu = cmplx(q_u(M_prime),     p_u(M_prime),     dp)
          du = cmplx(q_u(M_prime+Z),   p_u(M_prime+Z),   dp)
          ud = cmplx(q_d(M_prime),     p_d(M_prime),     dp)
          dd = cmplx(q_d(M_prime+Z),   p_d(M_prime+Z),   dp)

          s1 = s1 + conjg(uu) * dd
          s2 = s2 + conjg(ud) * du

          n_uu = n_uu + real(uu*conjg(uu), dp)
          n_ud = n_ud + real(ud*conjg(ud), dp)
          n_du = n_du + real(du*conjg(du), dp)
          n_dd = n_dd + real(dd*conjg(dd), dp)
        end do

        Rxx(i) = Rxx(i) + real(s1 + s2, dp)
        Rxy(i) = Rxy(i) + aimag(s1 + s2)
        Rzz(i) = Rzz(i) + (n_uu - n_ud - n_du + n_dd)
      end do
    end do
    !$OMP END DO


    deallocate(q_u, p_u, q_d, p_d, qtmp, ptmp, Hq, Hp)
  end block
  !$OMP END PARALLEL

  Rxx = Rxx / real(D,   dp)
  Rxy = Rxy / real(D,   dp)
  Rzz = Rzz / real(2*D, dp)

  t1 = omp_get_wtime()

  call write_two_col('QM_Rxx_n16_fast.dat', t, Rxx, NSTEPS)
  call write_two_col('QM_Rxy_n16_fast.dat', t, Rxy, NSTEPS)
  call write_two_col('QM_Rzz_n16_fast.dat', t, Rzz, NSTEPS)
  
  print '(A,F8.3,A)', 'Done. Elapsed ', t1-t0, ' s'
  print *, 'Output files:'
  print *, '  QM_Rxx_n16_fast.dat'
  print *, '  QM_Rxy_n16_fast.dat'
  print *, '  QM_Rzz_n16_fast.dat'

contains
  pure subroutine init_ket(M, is_up, q, p)
    integer,  intent(in)  :: M
    logical,  intent(in)  :: is_up
    real(dp), intent(out) :: q(:), p(:)
    integer :: idx
    q = 0.0_dp;  p = 0.0_dp
    if (is_up) then
      idx = M
    else
      idx = M + Z
    end if
    q(idx) = 1.0_dp
  end subroutine init_ket

  pure subroutine propag_4stage(q, p, dt, acoef, bcoef, a, omega, qtmp, ptmp, Hq, Hp)
    real(dp), intent(inout) :: q(:), p(:)
    real(dp), intent(in)    :: dt, acoef(:), bcoef(:), a(:), omega
    real(dp), intent(inout) :: qtmp(:), ptmp(:), Hq(:), Hp(:)
    integer :: j, m
    m = size(acoef)
    qtmp = q;  ptmp = p
    do j = 1, m
      call H_apply(qtmp, Hq, a, omega)
      ptmp = ptmp - bcoef(j)*dt*Hq
      call H_apply(ptmp, Hp, a, omega)
      qtmp = qtmp + acoef(j)*dt*Hp
    end do
    q = qtmp;  p = ptmp
  end subroutine propag_4stage

  pure subroutine H_apply(state, result, a, omega)
    real(dp), intent(in)  :: state(:)
    real(dp), intent(out) :: result(:)
    real(dp), intent(in)  :: a(:), omega
    integer(int64) :: st, partner, ebit, kbit
    integer :: k
    real(dp) :: se, sm, diag, amp
    integer :: nD
    nD = size(state)
    result = 0.0_dp
    ebit   = ishft(1_int64, n_spins)
    do st = 0_int64, int(nD-1, int64)
      amp = state(st+1)
      if (abs(amp) < 1.0e-12_dp) cycle
      se = merge(+1.0_dp, -1.0_dp, .not. btest(st, n_spins))
      diag = se * (omega*0.5_dp)
      do k = 1, n_spins
        if (btest(st, n_spins) .neqv. btest(st, k-1)) then
          kbit    = ishft(1_int64, k-1)
          partner = ieor(st, ior(ebit, kbit))
          result(partner+1) = result(partner+1) + 0.5_dp * a(k) * amp
        end if
        sm   = merge(+1.0_dp, -1.0_dp, .not. btest(st, k-1))
        diag = diag + se * 0.25_dp * a(k) * sm
      end do
      result(st+1) = result(st+1) + diag * amp
    end do
  end subroutine H_apply

  subroutine init_hyperfine(a)
    real(dp), intent(out) :: a(:)

    if (size(a) == 16) then
      a = [ -0.999985_dp, -0.7369246_dp, 0.511210_dp, -0.0826998_dp, &
            0.0655341_dp, -0.562082_dp, -0.905911_dp, 0.357729_dp, &
            0.358593_dp, 0.869386_dp, -0.232996_dp, 0.0388327_dp, &
            0.661931_dp, -0.930856_dp, -0.893077_dp, 0.0594001_dp ]
    else if (size(a) == 4) then
      a = [ -0.999985_dp, -0.7369246_dp, 0.511210_dp, -0.0826998_dp ]
    else
      stop 'not supported.'
    end if
  end subroutine init_hyperfine

  subroutine write_two_col(fname, x, y, n)
    character(*), intent(in) :: fname
    real(dp),     intent(in) :: x(:), y(:)
    integer,      intent(in) :: n
    integer :: u,i
    open(newunit=u, file=fname, status='replace', action='write')
    write(u,'(A)') '# Time    Value'
    do i=1,n
      write(u,'(2ES16.8)') x(i), y(i)
    end do
    close(u)
  end subroutine write_two_col
end program QM_nN_fast
