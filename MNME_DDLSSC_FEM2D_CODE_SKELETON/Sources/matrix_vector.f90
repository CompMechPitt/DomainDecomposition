module matrix_vector
  use matrix_class
  use vector_class
  use mpi
#ifdef ENABLE_OPENMP
	use omp_lib 
#endif
  
  implicit none

contains

  subroutine matvec ( A, x, y )
    implicit none
    ! Parameters 
    type(matrix), intent(in)    :: A 
    type(vector), intent(in)    :: x
    type(vector), intent(inout) :: y

    ! Locals
    integer              :: ierr,i,j

    ! *** Master on Numerical Methods in Engineering (MNME) ***
    ! *** Module Name: Domain Decomposition and Large Scale Scientific Computing (DDLSSC) ***
    ! *** TASK #5
#ifdef ENABLE_OPENMP
	real(8), allocatable :: temp(:)
#endif

	y%data = 0.0

#ifdef ENABLE_OPENMP
	!$omp parallel default(none) shared(A,x,y) private(I,J,temp)
	allocate(temp(A%m))
	temp = 0.0
	!$omp do schedule(static)
	! Each thread accesses N/P columns of A and N/P elements of x
	do j=1, A%n
		do i=1, A%m
			temp(i) = temp(i) + A%data(i,j)*x%data(j)
		end do
	end do
	!$omp end do
	!$omp critical
    y%data(:) = y%data(:) + temp(:)
	!$omp end critical
	
	deallocate(temp)
	!$omp end parallel 
	
		
#else
	do j=1, A%n
		do i=1, A%m
			y%data(i) = y%data(i) + A%data(i,j)*x%data(j)
		end do
	end do
#endif
  
  end subroutine matvec

end module matrix_vector
