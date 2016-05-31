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
    integer              :: ierr, I, J
    
#ifdef ENABLE_OPENMP
	real(8), allocatable :: aux(:)
#endif

	y%data = 0.0

#ifdef ENABLE_OPENMP
	!$omp parallel default(none) shared(A,x,y) private(I,J,aux)
	allocate(aux(A%m))
	aux = 0.0
	!$omp do schedule(static)
	! Each thread accesses N/P columns of A and N/P elements of x
	do J=1, A%n
		do I=1, A%m
			aux(I) = aux(I) + A%data(I,J)*x%data(J)
		end do
	end do
	!$omp end do
	!$omp critical
    y%data = y%data + aux
	!$omp end critical
	deallocate(aux)
	!$omp end parallel 
#else
	do J=1, A%n
		do I=1, A%m
			y%data(I) = y%data(I) + A%data(I,J)*x%data(J)
		end do
	end do
#endif

  end subroutine matvec

end module matrix_vector
