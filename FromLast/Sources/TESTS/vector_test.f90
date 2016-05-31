program vector_test
  use mpi
  implicit none

  type vector
     integer                  :: n         ! Number of vector entries
     real(8), allocatable     :: data(:)   ! Vector entries
     !type(data_desc), pointer :: desc      ! Pointer to the FE mesh descriptor
  end type vector
  
!  type data_desc
!    integer  :: nx_g, ny_g 
!    integer  :: nx_l, ny_l
!    integer  :: node_num
!    integer  :: elem_num

!    integer  :: me, np 

!    integer  :: up_pid
!    integer  :: up_off   ! Node identifier of the first node shared by this process
!                           ! and the one inmediately on top of it

!    integer  :: down_pid     
!    integer  :: down_off ! Node identifier of the first node shared by this process
!                           ! and the one inmediately below it
      
!    integer  :: mpi_comm
! end type data_desc
 
 type (vector) x
 type (vector) y
 real(8) :: alpha, alpha_f
 integer :: i, me, np, ierr, end

! type (data_desc) desc
 
 call MPI_Init(ierr)
 
 call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD, me, ierr)
 
! x%desc%me=me
! x%desc%np=np
! x%desc%mpi_comm=MPI_COMM_WORLD
 
! y%desc%me=me
! y%desc%np=np
! y%desc%mpi_comm=MPI_COMM_WORLD

 x%n=5
 y%n=5
 
 ALLOCATE(x%data(x%n),y%data(y%n))
 
 x%data=(/ 1, 1, 2, 1, 1/)
 y%data=(/ 2, 2, 3, 2, 2/)
 
! print *,me,'El vector x es:'
! DO i=1,x%n
!	print *,me,' ',x%data(i)
! END DO
 
!  print *,me,'El vector y es:'
! DO i=1,y%n
!	print *,me,' ',y%data(i)
! END DO
 
 !call vector_dot(x,y,alpha)
	
	do i=1, x%n
	
		alpha=alpha+x%data(i)*y%data(i)
	
	end do

!    if (me == 0) then
    
!		call MPI_Reduce(MPI_IN_PLACE,alpha,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		
!	else
	
!		call MPI_Reduce(alpha,MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     
!   end if
    
		
	print *,me,' antes ',alpha
		
		call MPI_Allreduce(alpha,alpha_f,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
			
	print *,me,' despues ',alpha_f
    
!    call mpi_barrier(x%desc%mpi_comm, ierr)
!    call mpi_abort(x%desc%mpi_comm, -1, ierr)
    
    call MPI_Finalize(ierr)

	
	print *,''
	! Impresión del vector por columnas
    !    DO i=1,x%n
    !            print *,alpha%data(i)
    !    END DO
    print *,me,'',alpha_f
	print *,''
	
	DEALLOCATE(x%data,y%data)

end program

!subroutine vector_dot ( x, y, alpha )
!    implicit none
!    ! Parameters 
!    type(vector), intent(in)  :: x
!    type(vector), intent(in)  :: y
!    real(8)     , intent(out) :: alpha
    
!    ! Locals
!    integer :: ierr, end, i

!    ! *** Master on Numerical Methods in Engineering (MNME) ***
!    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
!    ! *** TASK #1
!    ! *** CODE DEVELOPMENT REQUIRED ***
!    !     THE BODY OF THE FOLLOWING 
!    !     SUBROUTINE MUST BE DEVELOPED
!    !     FOR THE FINAL ASSIGNMENT
    
!	if (x%desc%me==x%desc%np-1) then
	
!		end = x%n/x%desc%np + mod(x%n,x%desc%np) 
	
!	else
	
!		end = x%n/x%desc%np
	
!	endif
	
!	do i=1, end
	
!		alpha=x%data(i)*y%data(i)
	
!	end do

!    if (x%desc%me == 0) then
    
!		call MPI_Reduce(MPI_IN_PLACE,alpha,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		
!	else
	
!		call MPI_Reduce(alpha,MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
     
!    end if
    
!    call mpi_barrier(x%desc%mpi_comm, ierr)
!    call mpi_abort(x%desc%mpi_comm, -1, ierr)

!  end subroutine vector_dot
