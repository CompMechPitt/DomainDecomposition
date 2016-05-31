module vector_class
  use mpi
  use data_desc_class
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector


  integer, parameter :: MAX_TAG      = 10000
  integer            :: current_tag  = 0 

contains
  
  subroutine vector_create ( desc, v )
    implicit none
    ! Parameters
    type(data_desc), intent(in), target :: desc
    type (vector)  , intent(out)        :: v

    ! Locals
    integer :: k, j, i 

    ! Point to the descriptor of the distributed FE mesh
    v%desc => desc 
    
    v%n = v%desc%node_num
    allocate ( v%data(1:v%n) )
    
    v%data = 0.0

  end subroutine vector_create

  subroutine vector_axpby ( alpha, beta, x, y )
!*****************************************************************************
!
!  vector_axpby performs the vector update y = beta*y + alpha*x,
!  where beta and alpha are scalars; x and y are vectors.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) alpha, beta, scalars involved in the
!    vector operator
!   
!    Input, type(vector) x. Input vector x.
!
!    Input/Output, type(vector) y. Updated vector y.
!
    implicit none
    ! Parameters
    real(8)     , intent(in)    :: alpha, beta
    type(vector), intent(in)    :: x
    type(vector), intent(inout) :: y

    y%data = beta * y%data + alpha * x%data
 
  end subroutine vector_axpby

  subroutine vector_copy ( x, y )
!*****************************************************************************
!
!  vector_copy copies x into y 
!
!  Parameters:
!
!    Input, type(vector) x. Input vector x.
!
!    Input/Output, type(vector) y. Vector y resulting from the copy.
!
    implicit none
    ! Parameters
    type(vector), intent(in)    :: x
    type(vector), intent(inout) :: y

    y%data = x%data
 
  end subroutine vector_copy

  subroutine vector_dot ( x, y, alpha )
    implicit none
    ! Parameters 
    type(vector), intent(in)  :: x
    type(vector), intent(in)  :: y
    real(8)     , intent(out) :: alpha
    ! real(8)     , intent(inout) :: temp_alpha
    ! Locals
    integer :: ierr
    integer :: i

    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #1
  
    
    do i=1, x%n
		alpha = alpha + x%data(i)*y%data(i)
	end do

	call MPI_Allreduce(mpi_in_place,&                 !send data
	alpha,                          &                 !receive data
	1,                              &                 !size    
	MPI_DOUBLE_PRECISION,           &                  !data type
	MPI_SUM,                        &                  !operation
	x%desc%mpi_comm,                 &              !communicator
	ierr) 
 
 
  end subroutine vector_dot


  subroutine vector_nrm2 ( x, nrm2 )
    implicit none
    ! Parameters 
    type(vector), intent(in)  :: x
    real(8)     , intent(out) :: nrm2

    ! Locals
    integer              :: ierr
     type(vector)  :: temp
    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #3
  
    nrm2 = 0.0
    
    
    call vector_create(x%desc,temp)
	call vector_comm(x,temp)

    call vector_dot(x,temp,nrm2)
    
    nrm2=sqrt(nrm2)
    
    call vector_destroy (temp)
    
  end subroutine vector_nrm2

  subroutine vector_comm ( x, y )
  !*****************************************************************************
!
!  vector_comm sums-up all subdomain contributions (partially 
!  summed entries) to vector entries of x shared among several 
!  subdomains, i.e., vector entries corresponding to FE mesh nodes
!  laying on the interface, and stores the sum into y. In the case 
!  we are facing (1D FE mesh data distribution), a given node on the 
!  interface is shared among at most 2 different subdomains.
!
!  Parameters:
!
!    Input, type(vector) x. Input (partially summed) vector x. 
!
!    Input/Output, type(vector) y. Output (fully summed) vector y.  
!
    implicit none
    ! Parameters 
    type(vector), intent(in)     :: x
    type(vector), intent(inout)  :: y
    
    ! Locals
    integer              :: status(MPI_STATUS_SIZE),ierr
    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #2
    
    !create variable to save the buffer for sending and receiving
    real(8), allocatable::recvtop(:)
    real(8), allocatable::recvbot(:)
    real(8), allocatable::sendtop(:)
    real(8), allocatable::sendbot(:)
  
   !allocate memory for buffer
    allocate(recvtop(1:x%desc%nx_l))
    allocate(recvbot(1:x%desc%nx_l))
    allocate(sendtop(1:x%desc%nx_l))
    allocate(sendbot(1:x%desc%nx_l))
    
     call vector_create(x%desc, y)
    sendtop = x%data(x%desc%up_off:x%desc%node_num)
    sendbot = x%data(1:x%desc%nx_l)
    
     call vector_copy(x,y)
     
     if (x%desc%np.gt.1) then
     !bottom to top 
		call MPI_Sendrecv(sendtop,  & ! send buffer
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION,	& ! data type
			x%desc%up_pid,			& ! destination ...
			99,				        & ! stag
			recvbot, 		    	& ! recv buffer
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION, 	& ! data type		
			x%desc%down_pid,  		& ! source ...
			99,  			      	& !  tag
			x%desc%mpi_comm,		& ! communicator
			status, ierr)
        !top to bottom
		call MPI_Sendrecv(sendbot,  & ! send buffer
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION,   & ! data type
			x%desc%down_pid,		& ! destination.
			99,			         	& ! stag
			recvtop, 			    & ! recv buffer
			x%desc%nx_l, 			& ! dimension
			MPI_DOUBLE_PRECISION,   & ! data type		
			x%desc%up_pid,  		& ! source ...
			99,  				    & !  tag
			x%desc%mpi_comm,		& ! communicator
			status, ierr)
			
			
	   if (x%desc%me == 0) then
            y%data(x%desc%up_off:x%desc%node_num) = y%data(x%desc%up_off:x%desc%node_num) + recvtop

        else if (x%desc%me == x%desc%np-1) then
            y%data(1:x%desc%nx_l) = y%data(1:x%desc%nx_l) + recvbot
	        
		else
            y%data(1:x%desc%nx_l) = y%data(1:x%desc%nx_l) + recvbot
            y%data(x%desc%up_off:x%desc%node_num) = y%data(x%desc%up_off:x%desc%node_num) + recvtop
        end if
    else
    
    write(*,*) 'Total number of process should be larger than 1,wrong in vector_comm'
    call mpi_barrier(x%desc%mpi_comm, ierr)
    call mpi_abort(x%desc%mpi_comm, -1, ierr)

    end if

    
    deallocate(recvtop)
    deallocate(recvbot)
    deallocate(sendtop)
    deallocate(sendbot)
   

  end subroutine vector_comm


  subroutine vector_weight ( x )
!*****************************************************************************
!
!  vector_weight multiplies each component of x shared among
!  several subdomains (i.e., vector entries corresponding to FE 
!  mesh nodes laying on the interface) by the inverse of the 
!  number of subdomains (i.e., by 1/2) 
!
!  Parameters:
!
!    Input/Output, type(vector) x. Output vector x.  
!
    implicit none
    ! Parameters 
    type(vector), intent(inout)  :: x 
    
    ! Locals
    integer              :: ierr

    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #4
  
    
    if (x%desc%np.gt.1) then
	
		if (x%desc%me /= x%desc%np-1) then
			x%data(x%desc%up_off:x%desc%node_num) = 0.5*x%data(x%desc%up_off:x%desc%node_num)
		endif

		if (x%desc%me /= 0) then
			x%data(1:x%desc%nx_l) = 0.5*x%data(1:x%desc%nx_l)
		endif
    else
    
    write(*,*) 'Total number of process should be larger than 1,wrong in vector_weight'
    call mpi_barrier(x%desc%mpi_comm, ierr)
    call mpi_abort(x%desc%mpi_comm, -1, ierr)

    end if
    
    

  end subroutine vector_weight


  subroutine vector_destroy ( v ) 
    implicit none
    ! Parameters
    type (vector), intent(inout)     :: v
   
    v%n = -1
    deallocate ( v%data )

  end subroutine vector_destroy

end module vector_class
