program vector_test2
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector
 
 type (vector) x
 real(8) :: alpha
 !integer, allocatable :: requests (:)
 !integer, dimension(MPI_STATUS_SIZE) :: status
 !integer, dimension(1) :: requests
 integer :: i, me, np, ierr, end, nx_l, up_off, node_num, TEST

  x%n=15
  nx_l=5
  up_off=11
  node_num=15
  
  allocate(x%data(x%n))
  
  x%data=1.0
 
  call MPI_Init(ierr)
 
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
  
  	if (me /= np - 1) then
		x%data(up_off:node_num) = 0.5*x%data(up_off:node_num)
	endif

	if (me /= 0) then
		x%data(1:nx_l) = 0.5*x%data(1:nx_l)
	endif

   call mpi_barrier(MPI_COMM_WORLD, ierr)
   call mpi_abort(MPI_COMM_WORLD, -1, ierr)

   call MPI_Finalize(ierr)

  deallocate(x%data)

end program
