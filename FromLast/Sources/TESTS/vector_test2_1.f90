program vector_test2
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector
 
 type (vector) y
 real(8) :: alpha
 real(8), allocatable :: aux(:), send(:), send2(:), send3(:)
 !integer, allocatable :: requests (:)
 integer, dimension(MPI_STATUS_SIZE) :: status
 !integer, dimension(1) :: requests
 integer :: i, me, np, ierr, end, nx_l, up_off, node_num, TEST, up_pid, down_pid

  y%n=15
  nx_l=5
  up_off=11
  node_num=15
  
  ALLOCATE(y%data(y%n),aux(nx_l),send(nx_l),send2(nx_l),send3(nx_l))
  !ALLOCATE(requests(np-1))
 
  call MPI_Init(ierr)
 
  call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
  
  
 
  y%data=(/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1/)
  
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
	
	
	
	!print *,'CUELGUE 2'
	
	!    print *,me
	!	DO i=1,nx_l
	!		print *,aux(i)
	!	END DO

	
	!call MPI_Ibarrier(MPI_COMM_WORLD, )
	
	if (me /= 0) then
	
	!call MPI_Waitall(np-1,1,status,ierr)
	
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
	
	if (me== 3) then

		DO i=1,nx_l
			print *,y%data(i)
		END DO
		
	endif
	
	print *,'Cambio procesador'
	
		if (me== 0) then

		DO i=1,nx_l
			print *,y%data(i)
		END DO
		
	endif


   call MPI_Finalize(ierr)
   


DEALLOCATE(y%data,aux,send,send2)

end program
