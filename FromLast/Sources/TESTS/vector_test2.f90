program vector_test2
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector
 
 type (vector) x
 type (vector) y
 real(8) :: alpha
 real(8), allocatable :: aux(:), send(:), send2(:), send3(:)
 integer, dimension(MPI_STATUS_SIZE) :: status
 integer :: i, me, np, ierr, end, nx_l, up_off, node_num, TEST

  x%n=15
  y%n=15
  nx_l=5
  up_off=11
  node_num=15
  
  ALLOCATE(x%data(x%n),y%data(y%n),aux(nx_l),send(nx_l),send2(nx_l),send3(nx_l))
 
  call MPI_Init(ierr)
 
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
  
  !allocate(statuses(MPI_STATUS_SIZE,np-1))
 
  x%data=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
  y%data=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
 
! print *,'El vector x es:'
! DO i=1,x%n
!	print *,x%data(i)
! END DO
 
!  print *,'El vector y es:'
! DO i=1,y%n
!	print *,y%data(i)
! END DO
     
!    ! Locals

    
!    ! *** Master on Numerical Methods in Engineering (MNME) ***
!    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
!    ! *** TASK #2

	y=x
	!call vector_copy(x, y)
	
	!send = x%data(up_off:node_num)
	!send2 = y%data(1:nx_l)
	
!  print *,'El vector send es:'
! DO i=1,nx_l
!	print *,send(i)
! END DO

	if (me/=np-1) then
	
		call MPI_Sendrecv(x%data(up_off:node_num), nx_l, MPI_DOUBLE_PRECISION, 1, 0, aux, nx_l, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
		
		!call MPI_Send(y%data(1:nx_l), nx_l, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, ierr)
		
	end if
	
	if (me /= 0) then
	
		do i = 1, nx_l
			send2(i) = send2(i) + aux(i)
		end do
		
		call MPI_Sendrecv(y%data(1:nx_l), nx_l, MPI_DOUBLE_PRECISION, 0, 0, aux, nx_l, MPI_DOUBLE_PRECISION, 1, 0, MPI_COMM_WORLD, status, ierr)
		
		!call MPI_Recv(send3, nx_l, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, status, ierr)
		
	end if
	
	if (me /= np-1) then
	
		y%data(up_off:node_num) = aux
		
	endif

   !call mpi_barrier(x%desc%mpi_comm, ierr)
   !call mpi_abort(x%desc%mpi_comm, -1, ierr)

   call MPI_Finalize(ierr)

DEALLOCATE(x%data,y%data,aux,send,send2)

end program
