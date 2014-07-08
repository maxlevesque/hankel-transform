
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! the module precision_kind defines the precision
! this is done in order to have a transferable program

module precision_kind
    implicit none
    integer(kind=kind(1)) , parameter :: i2b = kind ( 1 ) ! simple precision integer 
    integer(kind=i2b)     , parameter :: i4b = 2_i2b * i2b ! simple precision integer
    integer(kind=i2b)     , parameter :: dp = kind ( 0.0d0 ) ! double precision real
    integer(kind=i2b)     , parameter :: sp = kind ( 0.0 ) ! simple precision real
end module precision_kind


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module numbers
    use precision_kind
    implicit none
    real(dp), parameter :: pi = acos ( -1.0_dp )
    real(dp), parameter :: twopi = 2.0_dp * pi
    real(dp), parameter :: twopisq = 2.0_dp * pi ** 2
    real(dp), parameter :: fourpi = 4.0_dp * pi
end module numbers


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program hankel

    use precision_kind ! defines precisions
    use numbers

    implicit none

    integer(i2b) :: nb_input_pts ! total number of line of input file dat.in
    real(dp) , allocatable , dimension(:) :: x_in , y_in ! input x and y
    real(dp) :: kr ! k*r
    integer(i2b) :: direction ! 1 for hankel transform, 2 for inverse hankel transform
    real(dp) :: x ! radial distance
    integer(i2b) :: i , j ! dummy
    real(dp) :: A ! area under the curve
    real(dp) :: yj0 , yj1 ! dummy
    real(dp) :: xmax ! maximum value for r or k
    integer(i2b) :: nstep ! number of points between 0 and xmax
    
    ! read the total number of input points in dat.in -> nb_input_pts
    call get_nb_input_pts ( nb_input_pts )
    
    ! allocate accordingly an array that will contains all x (x_in) and another that will contain all y (y_in)
    allocate ( x_in ( nb_input_pts ) )
    allocate ( y_in ( nb_input_pts ) )
    
    ! read input dat.in
    open ( unit = 10 , file = 'dat.in' )
    do i = 1 , nb_input_pts
    read ( 10 , * ) x_in ( i ) , y_in ( i )
    end do
    
    close ( 10 )
    
    ! ask user if it want a transform (r -> k) or invers transform (k -> r)
    77 write (*,*) 'do you want a transform (r->k) or inverse transform (k->r) ? (1/2)'
    read (*,*) direction
    
    ! only correct input values are 1 and 2
    if ( direction /= 1 .and. direction /=2 ) then
    write (*,*) 'you can only type 1 or 2'
    ! send user back to choice 1 or 2
    go to 77
    end if
    
    ! open output file transformed.out and write a first comment line
    open ( unit = 11 , file = 'transformed.out' )
    if ( direction == 1 ) then
    write ( 11 , * ) '# k   f(k)'
    else if ( direction == 2 ) then
    write ( 11 , * ) '# r   f(r)'
    end if
    
    ! ask user the range over which he wants to compute f(x) (xmax)
    99 if ( direction == 1 ) then
    write ( * , * ) 'what is the maximum value you want for | k | (type a real number) ?'
    else if ( direction == 2 ) then
    write ( * , * ) 'what is the maximum value you want for | r | (type a real number) ?'
    end if
    read ( * , * ) xmax
    
    ! check validy of xmax
    
    if ( xmax <= 0.0_dp ) then
    write ( * , * ) 'xmax should be positive'
    go to 99
    end if
    
    ! ask user how many points he wants between 0 and xmax
    100 write ( * , * ) 'how many points do you want between 0 and ' , xmax , ' ? (enter an positive integer) '
    read ( * , * ) nstep
    
    ! check validy of xmax
    if ( nstep <= 0 ) then
    write ( * , * ) 'nstep should be positive'
    go to 100
    end if
    
    ! for nb_input_pts, one has nb_input_pts - 1 intervals on which to compute the trapeze area under the curve (A)
    ! for each x in between 0 and rmax (around 20) with a step of rstep (0.01)
    ! sum all A to get the total area which is the integral (hankel_int)
    do i = 0 , nstep
        ! for now x ranges from 0 to 20
        x = real ( i , kind = dp ) * xmax / real ( nstep , kind = dp )
        ! init integral f(x) = int ... f(k) ... dk
        A = 0.0_dp
        ! if x is zero then A = 0.0_dp. don't waste your time in the loop at dividing by zero. 
        if ( abs ( x ) <= epsilon ( 1.0_dp ) ) go to 88
    
        do j = 1 , nb_input_pts - 1
            if ( abs ( x_in ( j ) ) <= epsilon ( 1.0_dp ) ) then! care of dividing by zero
                yj0 = y_in ( j ) * x_in ( j ) ** 2 ! because sin(0)/0 = 1.0
            else 
                yj0 = y_in ( j ) * x_in ( j ) ** 2 * sin ( x_in ( j ) * x ) / ( x_in ( j ) * x )
            end if
            if ( abs ( x_in ( j + 1 ) ) <= epsilon ( 1.0_dp ) ) then
                yj1 = y_in ( j + 1 ) * x_in ( j + 1 ) ** 2 ! because sin(0)/0 = 1.0
            else
                yj1 = y_in ( j + 1 ) * x_in ( j + 1 ) ** 2 * sin ( x_in ( j + 1 ) * x ) / ( x_in ( j + 1 ) * x )
            end if
            A = A + ( x_in ( j + 1 ) - x_in ( j ) ) * ( yj0 + yj1 ) / 2.0_dp ! trapeze method for integral calculation
        end do
    
        ! normalize depending of transformation direction
        if ( direction == 1 ) then 
            A = A * fourpi
        else if ( direction == 2 ) then
            A = A / twopisq
        end if
    
        ! write output file x , f(x)
        88 write (11,*) x , A
    
    end do
    
    ! close output file
    close ( 11 )
    
    ! deallocate x_in and y_in
    deallocate ( x_in , y_in )

end program hankel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine computes the total number of lines of input file

subroutine get_nb_input_pts ( nb_input_pts )

    use precision_kind
    
    implicit none
    
    integer(i2b) :: ios ! iostat of the read statement: 
    integer(i2b) :: nb_input_pts

!    If the value of IOstatus is zero, the previous READ was executed flawlessly and all variables have received their input values. This is the normal case.
!    If the value of IOstatus is positive, the previous READ has encountered some problem. In general, without knowing the system dependent information, it is impossible to determine what the problem was. However, if hardware and I/O devices are working, a commonly seen problem would be illegal data. For example, supplying a real number to an integer variable. If IOstatus is positive, you cannot trust the values of the variables in the READ statement; they could all contain garbage values, or some of them are fine while the others are garbage.
!    If the value of IOstatus is negative, it means the end of the input has reached. Under this circumstance, some or all of the variables in the READ may not receive input values. 

    open ( 11 , file = 'dat.in' )

    ! init the total number of input points
    nb_input_pts = 0
    do while ( .true. )
        read ( 11 , * , iostat = ios )
        if ( ios > 0 ) then
            write (*,*) 'Error in compute_ck_dipolar.f90'
            write (*,*) 'something went wrong. stop'
            stop
        else if ( ios < 0 ) then
            exit! end of file reached
        else
            nb_input_pts = nb_input_pts + 1! increment total number of input points
        end if
    end do
    close ( 11 )

    ! if 0 lines then the input file is missing
    if ( nb_input_pts == 0 ) stop 'critical stop : input file dat.in is missing'

    ! inform user
    print*,'dat.in has ' , nb_input_pts , ' lines'

end subroutine get_nb_input_pts
