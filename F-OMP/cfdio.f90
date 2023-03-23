module cfdio

  implicit none

contains

subroutine writedatafiles(psi, m, n, scale)

  integer :: m, n, scale
  double precision ::  psi(0:m+1, 0:n+1)

  double precision, allocatable :: vel(:,:,:)
  integer, allocatable :: rgb(:,:,:)

  double precision :: modvsq, hue
  integer :: i, j, k

  integer, parameter :: iounitvel = 10, iounitcol = 11

! Compute local velocities and colours

  allocate(rgb(3,m,n))
  allocate(vel(2,m,n))

  do i = 1, m
     do j = 1, n

        vel(1,i,j) =   (psi(i,j+1)-psi(i,j-1)) / 2.0
        vel(2,i,j) = - (psi(i+1,j)-psi(i-1,j)) / 2.0

        modvsq = vel(1,i,j)**2 + vel(2,i,j)**2
        hue = modvsq**0.4

        call hue2rgb(hue, rgb(1,i,j), rgb(2,i,j), rgb(3,i,j))

     end do
  end do

!  Write out

  open(unit=iounitcol, file='colourmap.dat', form='formatted')
  open(unit=iounitvel, file='velocity.dat',  form='formatted')

  do j = 1, n
     do i = 1, m

!  Write colour map of velocity magnitude at every point

        write(iounitcol,fmt='(i4,1x,i4,1x,i3,1x,i3,1x,i3)') &
              i, j, rgb(1,i,j), rgb(2,i,j), rgb(3,i,j)

!  Only write velocity vectors every "scale" points
           
        if (mod(i-1,scale) == (scale-1)/2 .and. &
            mod(j-1,scale) == (scale-1)/2         ) then

           write(iounitvel,fmt='(i4,1x,i4,1x,g12.5,1x,g12.5)') &
                 i, j, vel(1,i,j), vel(2,i,j)
        end if
        
     end do
  end do

  close(unit=iounitcol)
  close(unit=iounitvel)

end subroutine writedatafiles


subroutine writeplotfile(m, n, scale)

  integer :: m, n, scale
  integer, parameter :: iounit = 10

  open(unit=iounit, file='cfd.plt', form='formatted')

  write(iounit,*) 'set size square'
  write(iounit,*) 'set key off'
  write(iounit,*) 'unset xtics'
  write(iounit,*) 'unset ytics'

  write(iounit,fmt='('' set xrange ['',i4,'':'',i4, '']'')') 1-scale, m+scale
  write(iounit,fmt='('' set yrange ['',i4,'':'',i4, '']'')') 1-scale, n+scale

  write(iounit,fmt='('' plot "colourmap.dat" w rgbimage, "velocity.dat" u 1:2:&
       &('',i2,''*0.75*$3/sqrt($3**2+$4**2)):&
       &('',i2,''*0.75*$4/sqrt($3**2+$4**2)) &
       &with vectors  lc rgb "#7F7F7F"'')') scale, scale

  close(unit=iounit)

end subroutine writeplotfile


subroutine hue2rgb(hue, r, g, b)

  double precision :: hue

  integer :: r, g, b
  integer, parameter :: rgbmax = 255

  r = rgbmax*colfunc(hue-1.0)
  g = rgbmax*colfunc(hue-0.5)
  b = rgbmax*colfunc(hue    )

end subroutine hue2rgb


double precision function colfunc(x)

  double precision :: x, absx, val

  double precision, parameter :: x1 = 0.2, x2 = 0.5

  absx = abs(x)

  if (absx .gt. x2) then
     val = 0.0
  else if (absx .lt. x1) then
     val = 1.0
  else
     val = 1.0 - ((absx-x1)/(x2-x1))**2
  end if

  colfunc = val
      
end function colfunc

double precision function gettime()

  logical, save :: firstcall = .true.

  integer, parameter :: int32kind = selected_int_kind( 9)
  integer, parameter :: int64kind = selected_int_kind(18)

  integer, parameter :: intkind = int64kind

  integer(kind = intkind) :: count,rate

  double precision, save :: ticktime

  if (firstcall) then

     firstcall = .false.

     call system_clock(count, rate)

     ticktime = 1.0d0/dble(rate)
     gettime  = dble(count)*ticktime

!     write(*,*) 'Clock resolution is ', ticktime*1.0e6, ', usecs'

  else

     call system_clock(count)

     gettime = dble(count)*ticktime

  end if

end function gettime

end module cfdio
