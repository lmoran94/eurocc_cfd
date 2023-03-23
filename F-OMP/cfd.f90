program cfd

  use boundary
  use jacobi
  use cfdio

  use omp_lib

  implicit none

! Output frequency
  
  integer, parameter :: printfreq = 1000

! Variables associated with convergence

  double precision :: error, bnorm

! Set tolerance for convergence; zero or negative means do not check

  double precision, parameter :: tolerance = 0.0d0

! Main arrays

  double precision, allocatable ::  psi(:,:), zet(:,:)
  double precision, allocatable ::  psitmp(:,:), zettmp(:,:)

! Command-line arguments

  integer :: scalefactor,  numiter

  double precision :: re  ! re = 3.7 seems to be stability limit with Jacobi

  integer, parameter :: maxline = 32
  character(len=maxline) :: tmparg

!  Basic sizes of simulation

  integer, parameter :: bbase = 10
  integer, parameter :: hbase = 15
  integer, parameter :: wbase =  5
  integer, parameter :: mbase = 32
  integer, parameter :: nbase = 32

  logical :: irrotational = .true., checkerr = .false.

!  Some auxiliary parameters and variables

  integer :: m, n, b, h, w
  integer :: iter

  integer :: nthread

  double precision :: tstart, tstop, ttot, titer, modvsq, hue

!  Are we stopping based on tolerance?

  if (tolerance .gt. 0.0) checkerr = .true.

!  Parallel initialisation

!  Read in parameters

  if (command_argument_count() /= 2 .and. command_argument_count() /= 3) then

     write(*,*) 'Usage: cfd <scale> <numiter> [reynolds]'
     stop

  end if

  call get_command_argument(1, tmparg)
  read(tmparg,*) scalefactor
  call get_command_argument(2, tmparg)
  read(tmparg,*) numiter

  if (command_argument_count() == 3) then

     irrotational = .false.
     call get_command_argument(3, tmparg)
     read(tmparg,*) re
        
  else

     re = -1.0
     
  end if

  if (.not. checkerr) then
     write(*,fmt='('' Scale factor = '',i3,'', iterations = '', i6)') &
           scalefactor, numiter
  else
     write(*,fmt='('' Scale factor = '',i3,'', iterations = '', i6, &
          &'', tolerance = '', g11.4)') scalefactor, numiter, tolerance
  end if

  if (irrotational) then
        
     write(*,*) 'Irrotational flow'
        
  else

     write(*,fmt='('' Reynolds number = '', f6.3)') re
        
  end if

!  Calculate b, h & w and m & n
        
  b = bbase*scalefactor 
  h = hbase*scalefactor
  w = wbase*scalefactor 
  m = mbase*scalefactor
  n = nbase*scalefactor

  re = re / dble(scalefactor)

  nthread = omp_get_max_threads()

  write(*,fmt='('' Running CFD on '', i4, '' x '', i4, '' grid using '', &
                  &i4, '' thread(s)'')') m, n, nthread

!  Allocate arrays, including halos on psi and tmp

  allocate(psi(0:m+1, 0:n+1))
  allocate(zet(0:m+1, 0:n+1))

  allocate(psitmp(0:m+1, 0:n+1))

  if (.not. irrotational) then

     allocate(zettmp(0:m+1, 0:n+1))

  end if

!  Zero the psi array

  psi(:,:) = 0.0
  zet(:,:) = 0.0

!  Set the psi boundary condtions which are constant

   call boundarypsi(psi, m, n, b, h, w)

!  Compute normalisation factor for error

   bnorm = sum(psi(:,:)**2)

   if (.not. irrotational) then

!    Update the zeta boundary condtions which depend on psi

     call boundaryzet(zet, psi, m, n)

!    Update the normalisation

     bnorm = bnorm + sum(zet(:,:)**2)

  end if

   bnorm = sqrt(bnorm)

!  Begin iterative Jacobi loop

   write(*,*)
   write(*,*) 'Starting main loop ...'
   write(*,*)

   tstart = gettime()

  do iter = 1, numiter

!  Compute the new psi based on the old one

     if (irrotational) then

!  Call function with no vorticity

        call jacobistep(psitmp, psi, m, n)

     else

!  Call function containing vorticity

        call jacobistepvort(zettmp, psitmp, zet, psi, m, n, re)

     end if

!  Compute current error value if required
     
     if (checkerr .or. iter == numiter) then

        error = deltasq(psitmp, psi, m, n)

        if (.not. irrotational) then

           error = error + deltasq(zettmp, zet, m, n)

        end if

        error = sqrt(error)
        
        error = error / bnorm

     end if

!  Copy back

!$omp parallel workshare
     psi(1:m, 1:n) = psitmp(1:m, 1:n)
!$omp end parallel workshare

     if (.not. irrotational) then

!$omp parallel workshare
        zet(1:m, 1:n) = zettmp(1:m, 1:n)
!$omp end parallel workshare

     end if

     if (.not. irrotational) then

!    Update the zeta boundary condtions which depend on psi

        call boundaryzet(zet, psi, m, n)
        
     end if

!  Quit early if we have reached required tolerance

     if (checkerr) then
        if (error .lt. tolerance) then
           write(*,*) 'CONVERGED iteration ', iter, ': terminating'
           exit
        end if
     end if

!  End iterative Jacobi loop

     if (mod(iter,printfreq) == 0) then

        if (.not. checkerr) then
           write(*,*) 'completed iteration ', iter
        else
           write(*,*) 'completed iteration ', iter, ', error = ', error
        end if

     end if

  end do

  if (iter .gt. numiter) iter = numiter

  tstop = gettime()

  ttot  = tstop-tstart
  titer = ttot/dble(iter)

  write(*,*) 
  write(*,*) '... finished'
  write(*,*)
  write(*,fmt='('' After    '', i6, '' iterations, error is '', g11.4)') &
        iter, error
  write(*,fmt='('' Time for '', i6, '' iterations was '',&
        &g11.4, '' seconds'')') iter, ttot
  write(*,fmt='('' Each individual iteration took '', g11.4, '' seconds'')') &
        titer
  write(*,*)
  write(*,*) 'Writing output file ...'

!  Output results

  call writedatafiles(psi, m, n, scalefactor)

!  Output gnuplot file

  call writeplotfile(m, n, scalefactor)

! Finish

  write(*,*) ' ... finished'
  write(*,*)
  write(*,*) 'CFD completed'
  write(*,*)

end program cfd

