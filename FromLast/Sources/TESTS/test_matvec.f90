  program test_matvec
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
  end type vector
  
  type matrix
     integer                  :: m         ! Number of rows in the matrix
     integer                  :: n         ! Number of cols in the matrix
     real(8), allocatable     :: data(:,:) ! Matrix entries
  end type matrix
  
 ! Locals
 integer              :: ierr, I, J, me, np, nx_l, up_off, node_num, up_pid, down_pid
 integer, dimension(MPI_STATUS_SIZE) :: status
 real(8), allocatable :: aux(:)

 ! Initialize matrix
 type (matrix) A
 A%m = 2
 A%n = 2
 A%data = 0.0
 A%data(1,1) = 1.0
 A%data(2,2) = 1.0
 
! Initialize fully-summed vector
 type (vector) x
 x%n = 2
 x%data=(/ 1, 1 /)

! Initialize aux vector and answer
 type (vector) y
 allocate(aux(A%m),y%data(A%m))
 aux = 0.0
 y%data = 0.0
 
! Initialize desc data
 nx_l = 1
 up_off = 2
 node_num = 2
 
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
   
	#ifdef ENABLE_OMP
		use omp_lib 
	#endif

	#ifdef ENABLE_OMP
	!$omp parallel default(none) shared(A,x,y) private(I,J)
	!$omp do reduction(+:y%data)
	! Each thread accesses N/P columns of A and N/P elements of x
	do J=1, A%n
		do I=1, A%m
			aux(I) = aux(I) + A%data(I,J)*x%data(J))
		end do
	end do
	!$omp end do
	!$omp critical comm
	!call vector_comm (aux,y)
	!$omp end critical comm
	!$omp end parallel 
	#else
	do J=1, A%n
		do I=1, A%m
			aux(I) = aux(I) + A%data(I,J)*x%data(J)
		end do
	end do
	!call vector_comm (aux,y)
	#endif
	
	print *,me,' ',aux(1),' ',aux(2)

    call mpi_barrier(x%desc%mpi_comm, ierr)
    call mpi_abort(A%desc%mpi_comm, -1, ierr)

  end program test_matvec
