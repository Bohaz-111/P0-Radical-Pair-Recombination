program QM_n4
  use iso_fortran_env
  implicit none

  integer, parameter :: dp = real64
  integer, parameter :: NSTEPS = 3000
  real(dp), parameter :: T_END = 30.0_dp
  real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp

  ! Hyperfine coupling constants (mT)
  real(dp), parameter :: a(4) = [-0.999985_dp, -0.7369246_dp, 0.511210_dp, -0.0826998_dp], omega = 0.5_dp

  integer, parameter :: D = 32
  integer, parameter :: Z = D

  real(dp), parameter :: integrator_a(4) = [0.5153528374311229_dp, -0.08578201941297365_dp, 0.4415830236164665_dp, 0.1288461583653842_dp]
  real(dp), parameter :: integrator_b(4) = [0.1344961992774311_dp, -0.2248198030794208_dp, 0.7563200005156683_dp, 0.3340036032863214_dp]

  complex(dp) :: H(D,D)

  complex(dp) :: q_vec(D), p_vec(D)
  complex(dp) :: q_temp(D), p_temp(D)

  complex(dp) :: U_up_up, U_up_down, U_down_up, U_down_down

  real(dp) :: t(NSTEPS), Rxx(NSTEPS), Rxy(NSTEPS), Rzz(NSTEPS)
  real(dp) :: dt
  integer :: i, M, M_prime
  
  call build_hamiltonian(a, H)
  
  dt = T_END / real(NSTEPS - 1, dp)
  
  do i = 1, NSTEPS
    t(i) = real(i - 1, dp) * dt
  end do
  
  Rxx = 0.0_dp
  Rxy = 0.0_dp
  Rzz = 0.0_dp
  
  do M = 1, D/2  
    call initialize_state_vector(M, .true., q_vec, p_vec)  
    call initialize_state_vector(M, .false., q_temp, p_temp) 
    
    do i = 1, NSTEPS
      call propagate_one_step(H, q_vec, p_vec, dt, integrator_a, integrator_b)
      call propagate_one_step(H, q_temp, p_temp, dt, integrator_a, integrator_b)
      
      do M_prime = 1, D/2
        call extract_matrix_elements(q_vec, p_vec, M_prime, .true., U_up_up, U_up_down)
        call extract_matrix_elements(q_temp, p_temp, M_prime, .false., U_down_up, U_down_down)
        
        call accumulate_correlations_single_time(U_up_up, U_up_down, U_down_up, U_down_down, &
                                               Rxx(i), Rxy(i), Rzz(i))
      end do
    end do
  end do
  
  Rxx = Rxx / real(Z, dp)
  Rxy = Rxy / real(Z, dp)
  Rzz = Rzz / real(2*Z, dp)

  call write_two_col('QM_Rxx_n4.dat', t, Rxx, NSTEPS)
  call write_two_col('QM_Rxy_n4.dat', t, Rxy, NSTEPS)
  call write_two_col('QM_Rzz_n4.dat', t, Rzz, NSTEPS)
  
  print *, 'Wrote QM_Rxx_n4.dat, QM_Rxy_n4.dat, QM_Rzz_n4.dat'




contains

  subroutine initialize_state_vector(M, is_up, q, p)
    integer, intent(in) :: M
    logical, intent(in) :: is_up
    complex(dp), intent(out) :: q(D), p(D)
    integer :: idx
    
    q = cmplx(0.0_dp, 0.0_dp, dp)
    p = cmplx(0.0_dp, 0.0_dp, dp)
    
    if (is_up) then
      idx = M  ! First half of basis states are up electron spins
    else
      idx = M + D/2  ! Second half of basis states are down electron spins
    end if
    
    q(idx) = cmplx(1.0_dp, 0.0_dp, dp)  ! Real part
  end subroutine initialize_state_vector

  subroutine propagate_one_step(H, q, p, dt, a_params, b_params)
    complex(dp), intent(in) :: H(D,D)
    complex(dp), intent(inout) :: q(D), p(D)
    real(dp), intent(in) :: dt
    real(dp), intent(in) :: a_params(4), b_params(4)
    
    complex(dp) :: q_temp(D), p_temp(D)
    complex(dp) :: Hq(D), Hp(D)
    integer :: j
    
    q_temp = q
    p_temp = p
    
    do j = 1, 4
      Hq = matmul(H, q_temp)
      p_temp = p_temp - b_params(j) * dt * Hq
      
      Hp = matmul(H, p_temp)
      q_temp = q_temp + a_params(j) * dt * Hp
    end do
    
    q = q_temp
    p = p_temp
  end subroutine propagate_one_step

  subroutine extract_matrix_elements(q, p, M_prime, is_up, U_up, U_down)
    complex(dp), intent(in) :: q(D), p(D)
    integer, intent(in) :: M_prime
    logical, intent(in) :: is_up
    complex(dp), intent(out) :: U_up, U_down
    
    integer :: idx_up, idx_down
    complex(dp) :: psi_up, psi_down
    
    idx_up = M_prime
    idx_down = M_prime + D/2
    
    psi_up = q(idx_up) + cmplx(0.0_dp, 1.0_dp, dp) * p(idx_up)
    psi_down = q(idx_down) + cmplx(0.0_dp, 1.0_dp, dp) * p(idx_down)
    U_up = psi_up
    U_down = psi_down

  end subroutine extract_matrix_elements

  subroutine accumulate_correlations_single_time(U_up_up, U_up_down, U_down_up, U_down_down, Rxx, Rxy, Rzz)
    complex(dp), intent(in) :: U_up_up, U_up_down, U_down_up, U_down_down
    real(dp), intent(inout) :: Rxx, Rxy, Rzz
    
    complex(dp) :: complex_correlation
    real(dp) :: real_part, imag_part
    
    complex_correlation = conjg(U_up_up) * U_down_down + conjg(U_up_down) * U_down_up
    real_part = real(complex_correlation, dp)
    imag_part = aimag(complex_correlation)
    
    Rxx = Rxx + real_part
    Rxy = Rxy + imag_part
    
    real_part = abs(U_up_up)**2 - abs(U_up_down)**2 - abs(U_down_up)**2 + abs(U_down_down)**2
    
    Rzz = Rzz + real_part
  end subroutine accumulate_correlations_single_time


  subroutine build_hamiltonian(a, H)
    real(dp), intent(in) :: a(4)
    complex(dp), intent(out) :: H(D,D)
    real(dp), parameter :: omega = 0.5_dp
    complex(dp) :: sx2(2,2), sy2(2,2), sz2(2,2), id2(2,2)
    complex(dp) :: tmp(4,4), tmp2(8,8), tmp3(16,16)
    complex(dp) :: Ixe(32,32), Iye(32,32), Ize(32,32)
    complex(dp) :: Sxloc(32,32), Syloc(32,32), Szloc(32,32)
    integer :: k
    
    id2 = reshape([(1.0_dp,0.0_dp), (0.0_dp,0.0_dp), (0.0_dp,0.0_dp), (1.0_dp,0.0_dp)], [2,2])
    sx2 = 0.5_dp * reshape([(0.0_dp,0.0_dp), (1.0_dp,0.0_dp), (1.0_dp,0.0_dp), (0.0_dp,0.0_dp)], [2,2])
    sy2 = 0.5_dp * reshape([(0.0_dp,0.0_dp), (0.0_dp,-1.0_dp), (0.0_dp,1.0_dp), (0.0_dp,0.0_dp)], [2,2])
    sz2 = 0.5_dp * reshape([(1.0_dp,0.0_dp), (0.0_dp,0.0_dp), (0.0_dp,0.0_dp), (-1.0_dp,0.0_dp)], [2,2])

    call build_operators(Sxloc, Syloc, Szloc)

    H = omega*Szloc
    
    do k = 1, 4
       select case(k)
       case(1)
          call kron2(id2, sx2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, id2, Ixe)
          call kron2(id2, sy2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, id2, Iye)
          call kron2(id2, sz2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, id2, Ize)
       case(2)
          call kron2(id2, id2, tmp);  call kron2(tmp, sx2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, id2, Ixe)
          call kron2(id2, id2, tmp);  call kron2(tmp, sy2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, id2, Iye)
          call kron2(id2, id2, tmp);  call kron2(tmp, sz2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, id2, Ize)
       case(3)
          call kron2(id2, id2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, sx2, tmp3); call kron2(tmp3, id2, Ixe)
          call kron2(id2, id2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, sy2, tmp3); call kron2(tmp3, id2, Iye)
          call kron2(id2, id2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, sz2, tmp3); call kron2(tmp3, id2, Ize)
       case(4)
          call kron2(id2, id2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, sx2, Ixe)
          call kron2(id2, id2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, sy2, Iye)
          call kron2(id2, id2, tmp);  call kron2(tmp, id2, tmp2); call kron2(tmp2, id2, tmp3); call kron2(tmp3, sz2, Ize)
       end select

       H = H + a(k) * (matmul(Sxloc, Ixe) + matmul(Syloc, Iye) + matmul(Szloc, Ize))
    end do
  end subroutine build_hamiltonian

  subroutine build_operators(Sx, Sy, Sz)
    complex(dp), intent(out) :: Sx(D,D), Sy(D,D), Sz(D,D)
    complex(dp) :: sx2(2,2), sy2(2,2), sz2(2,2), id2(2,2)
    complex(dp) :: tmp(4,4), tmp2(8,8), tmp3(16,16)

    id2 = reshape([(1.0_dp,0.0_dp), (0.0_dp,0.0_dp), (0.0_dp,0.0_dp), (1.0_dp,0.0_dp)], [2,2])
    sx2 = 0.5_dp * reshape([(0.0_dp,0.0_dp), (1.0_dp,0.0_dp), (1.0_dp,0.0_dp), (0.0_dp,0.0_dp)], [2,2])
    sy2 = 0.5_dp * reshape([(0.0_dp,0.0_dp), (0.0_dp,-1.0_dp), (0.0_dp,1.0_dp), (0.0_dp,0.0_dp)], [2,2])
    sz2 = 0.5_dp * reshape([(1.0_dp,0.0_dp), (0.0_dp,0.0_dp), (0.0_dp,0.0_dp), (-1.0_dp,0.0_dp)], [2,2])

    call kron2(sx2, id2, tmp)              ! 4×4
    call kron2(tmp, id2, tmp2)             ! 8×8
    call kron2(tmp2, id2, tmp3)            ! 16×16
    call kron2(tmp3, id2, Sx)              ! 32×32

    call kron2(sy2, id2, tmp);  call kron2(tmp, id2, tmp2);  call kron2(tmp2, id2, tmp3);  call kron2(tmp3, id2, Sy)
    call kron2(sz2, id2, tmp);  call kron2(tmp, id2, tmp2);  call kron2(tmp2, id2, tmp3);  call kron2(tmp3, id2, Sz)
  end subroutine build_operators

  subroutine kron2(A, B, C)
    complex(dp), intent(in) :: A(:,:), B(:,:)
    complex(dp), intent(out) :: C(size(A,1)*size(B,1), size(A,2)*size(B,2))
    integer :: i,j, p,q, m1,n1,m2,n2
    m1=size(A,1); n1=size(A,2); m2=size(B,1); n2=size(B,2)
    do i=1,m1
      do j=1,n1
        do p=1,m2
          do q=1,n2
            C((i-1)*m2+p,(j-1)*n2+q) = A(i,j) * B(p,q)
          end do
        end do
      end do
    end do
  end subroutine kron2

  subroutine write_two_col(fname, x, y, n)
    character(*), intent(in) :: fname
    real(dp), intent(in) :: x(n), y(n)
    integer, intent(in) :: n
    integer :: u,i
    open(newunit=u, file=fname, status='replace', action='write')
    do i=1,n
      write(u,'(2ES20.12)') x(i), y(i)
    end do
    close(u)
  end subroutine write_two_col

end program QM_n4
