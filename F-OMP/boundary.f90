module boundary

  implicit none

contains

subroutine boundarypsi(psi, m, n, b, h, w)

  integer :: m, n, b, h, w
  
  double precision, dimension(0:m+1, 0:n+1) :: psi

  integer :: i, j

!  Set the boundary conditions on the bottom edge

  do i = b+1, b+w-1
     psi(i, 0) = float(i-b)
  end do

  do i = b+w, m
     psi(i, 0) = float(w)
  end do

  !  Set the boundary conditions on the right hand side

  do j = 1, h

     psi(m+1,j) = float(w)
     
  end do

  do j = h+1, h+w-1

     psi(m+1,j) = float(w-j+h)

  end do

end subroutine boundarypsi

subroutine boundaryzet(zet, psi, m, n)

  integer :: m, n
  
  double precision, dimension(0:m+1, 0:n+1) :: zet, psi

  integer :: i, j

! Set the zeta boundary conditions which depend on psi

  do j = 1, n

     zet(0,  j) = 2.0*(psi(1,j) - psi(0,  j))
     zet(m+1,j) = 2.0*(psi(m,j) - psi(m+1,j))

  end do

  do i = 1, m
     zet(i,0) = 2.0*(psi(i,  1)-psi(i,0))
  end do

  do i = 1, m
     zet(i,n+1) = 2.0*(psi(i,n)-psi(i,n+1))
  end do

end subroutine boundaryzet

end module boundary
