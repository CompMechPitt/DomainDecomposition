program precond_test
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector
  
  type matrix
     integer                  :: m         ! Number of rows in the matrix
     integer                  :: n         ! Number of cols in the matrix
     real(8), allocatable     :: data(:,:) ! Matrix entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type matrix
  
  type precond
     !type(data_desc), pointer :: desc        ! Pointer to the FE mesh descriptor
     integer				  :: node_num
     integer                  :: type
     type(vector)             :: inv_diag    ! Stores the inverse of the diagonal of A
     real(8), allocatable     :: chol(:,:)   ! Stores the cholesky factor of A in its 
                                             ! upper triangle
  end type precond
  
    integer, parameter  :: identity = 0  ! M = I
	integer, parameter  :: diagonal = 1  ! M = D, where D=diag(A)
	integer, parameter  :: neumann  = 2  ! M = sum_i=1^P = I_i A_i^-1 I_i^T
  
    type (matrix) A
    type (precond) M

    ! Locals
    integer              :: ierr, k, j, i, info, np, me

	A%m = 3
	A%n = 3
	M%node_num = A%m
	
	allocate(A%data(A%m,A%n))
	allocate(M%inv_diag%data(M%node_num))
	
	A%data = reshape((/1,4,7,2,5,8,3,6,9/),shape(A%data))
	
	call MPI_Init(ierr)
 
	call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
	call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
	
	DO i=1,A%m
		DO j=1,A%n
			print *,A%data(i,j)
		END DO
    END DO


    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #6
    
    ! Point to the descriptor of the distributed FE mesh
!    M%desc => desc
    
    M%type = type ! M = I nothing else should be done 
	
	if type == diagonal ! M = D^{-1}, where D=diag(A)
	
!		call vector_create(desc,M%inv_diag)
		
		do i = 1,M%node_num
			M%inv_diag%data(i) = A%data(i,i)
		end do
		
		!call vector_comm(M%inv_diag,M%inv_diag)
		if (me /= 0 .and. me /= np-1) then
	
		up_pid = me + 1
		down_pid = me - 1
					
	elseif (me == np-1) then
	
		up_pid = 0
		down_pid = me -1		
	
	elseif (me == 0) then
	
		up_pid = me + 1
		down_pid = np-1
		
	endif
  
   
!  print *,'El vector y es:'
! DO i=1,y%n
!	print *,y%data(i)
! END DO


! DO i=1,nx_l
!	print *,send(i)
! END DO
		
    call MPI_Sendrecv(y%data(up_off:node_num), &
						  nx_l, MPI_DOUBLE_PRECISION, &
						  up_pid, 0, aux, nx_l, &
						  MPI_DOUBLE_PRECISION, up_pid, &
						  0, MPI_COMM_WORLD, status, ierr)
	
	print *,'CUELGUE'

	
	if (me /= 0) then
	
		do i = 1, nx_l
			y%data(i) = y%data(i) + aux(i)
		end do
		
	endif
		
		call MPI_Sendrecv(y%data(1:nx_l), nx_l, &
						  MPI_DOUBLE_PRECISION, down_pid, &
						  0, aux, nx_l, MPI_DOUBLE_PRECISION, &
						  down_pid, 0, MPI_COMM_WORLD, &
						  status, ierr)
		

	
	!call MPI_Waitall(np-1,requests,status,ierr)

	
	if (me /= np-1) then
	
		y%data(up_off:node_num) = aux
		
	endif
		
		
		
		do i = 1,M%node_num
			M%inv_diag%data(i) = 1/M%inv_diag%data(i)
		end do
	
!	elseif type == neumann ! M = sum_i=1^P = I_i A_i^-1 I_i^T
	
!		allocate(M%chol(A%m,A%n))
		
!		M%chol = A%data
		
!		call DPOTRF('U',A%n,M%chol,A%m,info)
	
	end


    !call mpi_barrier(A%desc%mpi_comm, ierr)
    !call mpi_abort(A%desc%mpi_comm, -1, ierr)
 
    
  call MPI_Finalize(ierr)   


DEALLOCATE(A%data)

end program

