program XYZ_exchange_tracker

    ! Declare variables
    implicit none
    double precision, parameter :: mH = 2.d0, mB = 10.8d0, mC = 12d0, mAz = 14.d0, mO = 16.d0, mK = 39.d0
    double precision, dimension(1:3) :: ABC, ABChalf
    double precision :: A, B, C, temp, mtot, temp_d2
    double precision, allocatable :: m_elem(:), r(:,:), r_noPBC(:,:), r_elem(:,:,:)
    integer, parameter :: max_Z_elem = 19
    integer :: natom, zatom, ielem, max_n_elem, istep, iatom, jatom, zjatom, i, j, iarg, k, temp_int, reason
    integer, allocatable :: n_elem(:), index_elem(:,:), instant_ielem(:), user(:)
    character(80) :: XYZ, SYSTEM, arg, temp_str
    character(500) :: comment_line
    character(3), allocatable :: element(:), elem_zatom(:)    
    !END Declare variables

    ! Variables description ! 
    ! zatom : atomic number
    ! iatom, natom : ith atom of the system that contains natom atoms
    ! ielem : ith atom of element zatom
    ! nH, nB... number of specified element in system
    ! elem_zatom(zatom) : name of element (e.g. "H")
    ! m_elem(zatom) : mass of element, double precision
    ! n_elem(zatom) : number of atoms of element zatom in system
    ! index_elem(zatom, ielem) : atom index over natom, for ielem atom of element zatom
    ! r_elem(zatom, ielem, k) : (x,y,z) coordinates of the ielem atom of element zatom
    ! instant_ielem(zatom) : temporary variable to follow ielem 
    ! END variables description
    
    ! WARNINGS
    print *, "NOTE: SIGSEV FAULT can be due to an atom type in your system not declared yet in the program. "
    print *, "To include a new element, modify max_Z_elem, m_elem, and elem_zatom"
    print *, " "
    ! END WARNINGS
  
    ! Read bash arguments
    call readArguments(A,B,C,natom,SYSTEM,XYZ)
    ABC=[A,B,C]
    ABChalf = ABC/2
    i = 0
    print *, 'SYSTEM = ', trim(SYSTEM)
    print *, 'XYZ = ', trim(XYZ)
    write(*, '(" ABC = ", f12.3, f12.3, f12.3)') ABC
    print *, 'natom = ', natom
    ! END Read bash arguments
    

    !write(*,'("ABC = ", f12.2)') A

    ! Setting up file reading !
    allocate( element(natom))
    allocate(r(natom,3), source = 666.d0)
    allocate(r_noPBC(natom,3), source = 666.d0)

    ! Reading step 1
    open(001, file = XYZ)
    read(001, *) !natom
    read(001, *) comment_line !comment line
    do iatom = 1, natom ! Getting elements and coordinates
        read(001, *) element(iatom), r(iatom, 1), r(iatom, 2), r(iatom, 3)
    end do
    ! END reading step 1

    ! Setting up elements names !
    allocate(elem_zatom(max_Z_elem), source = '   ')
    elem_zatom(1) = "H"
    elem_zatom(5) = "B"
    elem_zatom(6) = "C"
    elem_zatom(7) = "N"
    elem_zatom(8) = "O"
    elem_zatom(19) = "K"
    ! END setting up elements names !

    ! Counting atoms per element type !
    allocate(n_elem(max_Z_elem), source = 0)
    do zatom = 1, max_Z_elem 
        n_elem(zatom) = count(element == elem_zatom(zatom))
    end do
    max_n_elem = maxval(n_elem) ! maximum number of atoms for one element type

    print *, "Detected number of atoms per element: "
    do zatom = 1, max_Z_elem
        if (n_elem(zatom) /= 0) print *, n_elem(zatom), elem_zatom(zatom)
    end do
    ! END counting atoms per element type !

    ! Setting up elements masses !
    allocate(m_elem(max_Z_elem) , source = 0.d0)
    m_elem(1) = 2.d0 ! Deuterium instead of H 
    m_elem(5) = 10.8d0 ! B
    m_elem(6) = 12.d0 ! C
    m_elem(7) = 14.d0 ! N
    m_elem(8) = 16.d0 ! O
    m_elem(19) = 39.d0 ! K
    ! END setting up elements masses !
    mtot = sum(n_elem*m_elem)
    write(*, '(" Total mass of system mtot = ", f12.3)') mtot

    ! Filling in index_elem to keep track of atomic indexes per element !
    allocate(index_elem(19,max_n_elem), source = -1)
    allocate(instant_ielem(max_Z_elem), source = 0)
    do iatom = 1, natom
        zatom = findloc(elem_zatom, element(iatom), 1)
        instant_ielem(zatom) = instant_ielem(zatom) + 1
        index_elem(zatom, instant_ielem(zatom)) = iatom
    end do
    ! END filling in index_elem !

    close(001)
    ! END setting up file reading


    ! Read XYZ file ! 
    print *, " "
    print *, "Starting analysis..."
    open(001, file = XYZ)
    open(101, file = trim(SYSTEM) // '_exchange_track.xyz')
    istep = 0
    do ! step loop
        istep = istep + 1
      allocate( user(natom), source = 0 )
      allocate(r_elem(max_Z_elem, max_n_elem, 3), source = 0.d0)
      deallocate(instant_ielem)
      allocate(instant_ielem(max_Z_elem), source = 0)

      ! Run time info 
      if( mod(istep,5000)==0 ) print*,"READING step ",istep
      ! END run time info

      read(001, *, IOSTAT=reason) !natom
      if (reason <0) then
        print *, "EOF was reached, program exiting normally"
        print *, " "
        print *, "Output file is ", trim(SYSTEM) // '_exchange_track.xyz'
        close(001)
        close(101)
        exit
      endif

      read(001, '(A)') comment_line !time line
      
      ! Read and store atomic coordinates ! 
      do iatom = 1, natom ! Getting elements and coordinates
          read(001, *, IOSTAT=reason) temp_str, r(iatom, 1), r(iatom, 2), r(iatom, 3)

          do k = 1, 3
              temp = r(iatom, k)
              r_noPBC(iatom, k) = temp ! Copy coordinates without PBC wrapping for the output
          end do
      end do

      do k = 1, 3 
          r(:, k) = r(:, k) - ABC(k) * floor(r(:, k)/ABC(k)) ! PBC wrapping 
      end do

      do zatom = 1, max_Z_elem
        do ielem = 1, n_elem(zatom)
          iatom = index_elem(zatom, ielem)
          r_elem(zatom, ielem, :) = r(iatom, :)
        end do
      end do
      ! END read and store atomic coordinates !

      ! Find closest O to every H
      zatom = 1
      zjatom = 8
      do i = 1, n_elem(zatom)
        temp = 10000.d0
        temp_d2 = temp**2
        temp_int = -1 
        iatom = index_elem(zatom, i)
        do j = 1, n_elem(zjatom)
            jatom = index_elem(zjatom, j)
            temp_d2 = dist2(iatom, jatom)
            if (temp_d2 < temp) then
                temp_int = j
                temp = temp_d2
            end if
        end do
        temp_int = index_elem(zjatom, temp_int)
        user(temp_int) = user(temp_int) + 1
      end do
      !END find closest O to every H

      ! Writing output ! 
      write(101, *) natom
      write(101, *) trim(comment_line)
      do iatom = 1, natom
        write(101, *) element(iatom), r_noPBC(iatom, :), user(iatom)
      end do
      ! END writing output !
    
    deallocate( user, r_elem )
    end do !step loop!
    
    !END Read XYZ file!

    ! subroutines & functions
    contains
    subroutine readArguments(A,B,C,natom,SYSTEM,XYZ)
        double precision, intent(out) :: A, B, C
        integer, intent(out) :: natom
        character (len=*), intent(out) :: SYSTEM, XYZ
        call get_command_argument(1,arg,status=iarg)
        read(arg,*) A
        call get_command_argument(2,arg,status=iarg)
        read(arg,*) B
        call get_command_argument(3,arg,status=iarg)
        read(arg,*) C
        call get_command_argument(4,arg,status=iarg)
        read(arg,*) natom
        call get_command_argument(5,arg,status=iarg)
        read(arg,*) SYSTEM
        call get_command_argument(6,arg,status=iarg)
        read(arg,*) XYZ
    end subroutine

    function dist2(iatom, jatom) 
        double precision :: dist2
        integer, intent(in) :: iatom, jatom
        double precision, dimension(1:3) :: dd
        dist2 = 0.d0
        dd = [0, 0, 0]
        do k = 1, 3
            dd(k) = abs( r(iatom, k) - r(jatom, k) )
            if ( dd(k) > ABChalf(k) ) dd(k) = ABC(k) - dd(k)
            dist2 = dist2 + dd(k)**2
        end do
    end function 
    ! END subroutines & functions

end program XYZ_exchange_tracker

  
