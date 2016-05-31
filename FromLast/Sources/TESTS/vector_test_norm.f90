program vector_test_norm
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector
 
 type (vector) vec_x, vec_y
 real(8) :: alpha, alpha_f, nrm2
 real(8), allocatable :: aux(:), send(:), send2(:), send3(:)
 !integer, allocatable :: requests (:)
 integer, dimension(MPI_STATUS_SIZE) :: status
 !integer, dimension(1) :: requests
 integer :: i, me, np, ierr, end, nx_l, up_off, node_num, TEST, up_pid, down_pid
	
  vec_x%n=15	
  vec_y%n=15
  nx_l=5
  up_off=11
  node_num=15
  
  ALLOCATE(vec_y%data(vec_y%n),vec_x%data(vec_x%n),aux(nx_l),send(nx_l),send2(nx_l),send3(nx_l))
  !ALLOCATE(requests(np-1))
 
  call MPI_Init(ierr)
 
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
  
 
  vec_x%data=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
  vec_y%data=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)

		!call comm(x,y)
		
		
		
		!call vector_dot(vec_x,vec_y,alpha)
		
		do i=1, vec_x%n
	
			alpha = alpha + vec_x%data(i)*vec_y%data(i)
	
		end do

		call MPI_Allreduce(alpha,alpha_f,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
		
		print *,me,' alpha_f ',alpha_f
		
		nrm2=(alpha_f)**0.5
		
		print *,me,' norm ',nrm2

  call MPI_Finalize(ierr)   


DEALLOCATE(vec_x%data,vec_y%data,aux,send,send2)

end program

 
