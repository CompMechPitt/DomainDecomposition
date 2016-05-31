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

  subroutine vector_dot ( x, y, alpha_f )
    implicit none
    ! Parameters 
    type(vector), intent(in)  :: x
    type(vector), intent(in)  :: y
    real(8)		, intent(out) :: alpha_f
    
    ! Locals
    integer :: ierr, i
    real(8) :: alpha
    
    alpha = 0.0

    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #1
	
	do i=1, x%n
		alpha = alpha + x%data(i)*y%data(i)
	end do

	call MPI_Allreduce(alpha,alpha_f,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

  end subroutine vector_dot


  subroutine vector_nrm2 ( x, nrm2 )
    implicit none
    ! Parameters 
    type(vector), intent(in)  :: x
    real(8)     , intent(out) :: nrm2

    ! Locals
    integer              :: ierr
    type(vector)		 :: y
    
    nrm2 = 0.0

	call vector_comm(x,y)
	call vector_dot(x,y,nrm2)
	
	nrm2=(nrm2)**0.5

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
    integer              :: ierr, status(MPI_STATUS_SIZE), i, end
    real(8), allocatable :: aux(:)
    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #2

	call vector_create(x%desc, y)

	call vector_copy(x,y)

	allocate(aux(y%desc%nx_l))
	aux = 0.0
	
	if (y%desc%np > 1) then
		
		call MPI_Sendrecv(y%data(y%desc%up_off:y%desc%node_num), &
							  y%desc%nx_l, MPI_DOUBLE_PRECISION, &
					         y%desc%up_pid, 0, aux, y%desc%nx_l, &
					      MPI_DOUBLE_PRECISION, y%desc%down_pid, &
					            0, MPI_COMM_WORLD, status, ierr)
	
		if (y%desc%me /= 0) then
	
			do i = 1, y%desc%nx_l
				y%data(i) = y%data(i) + aux(i)
			end do
		
		endif
	
		call MPI_Sendrecv(y%data(1:y%desc%nx_l), y%desc%nx_l, &
					   MPI_DOUBLE_PRECISION, y%desc%down_pid, &
				   0, aux, y%desc%nx_l, MPI_DOUBLE_PRECISION, &
					        y%desc%up_pid, 0, MPI_COMM_WORLD, &
					                            status, ierr)
	
		if (y%desc%me /= y%desc%np-1) then
	
			y%data(y%desc%up_off:y%desc%node_num) = aux
		
		endif

	end if

	deallocate(aux)

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
