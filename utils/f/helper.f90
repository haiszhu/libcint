      ! compute Vijkl tensor from Vmunu and idcoefs
      subroutine computeVijkl(nd, Norb, idcoefs, Vmunu, Vijkl)
      implicit none
      integer, intent(in) :: nd, Norb
      real *8, intent(in) :: idcoefs(nd, Norb*(Norb+1)/2)
      real *8, intent(in) :: Vmunu(nd, nd)
      real *8, intent(inout) :: Vijkl(Norb, Norb, Norb, Norb)

      ! Local variables
      real *8, allocatable :: collocation_matrix(:, :)
      real *8, allocatable :: collocation_matrix_ij(:)
      real *8, allocatable :: collocation_matrix_kl(:)
      real *8, allocatable :: outer_product(:, :)
      integer :: i, j, k, l, tmpidx

      ! Allocate collocation_matrix (nd × Norb^2)
      allocate(collocation_matrix(nd, Norb**2))

      ! Construct collocation_matrix
      tmpidx = 0
      do i = 1,Norb
        do j = 1,i
          tmpidx = tmpidx + 1
          collocation_matrix(:,(j-1)*Norb+i) = idcoefs(:,tmpidx)
        enddo
      enddo

      ! Symmetrize collocation_matrix
      do i = 1,Norb
        do j = i+1,Norb
          collocation_matrix(:,(j-1)*Norb+i) = &
              collocation_matrix(:,(i-1)*Norb+j)
        enddo
      enddo

      ! Zero out Vijkl
      Vijkl = 0.0d0

      ! Allocate temporary vectors
      allocate(collocation_matrix_ij(nd), collocation_matrix_kl(nd))
      allocate(outer_product(nd, nd))

      ! Compute Vijkl
      !$omp parallel do &
      !$omp private(i, j, k, l, collocation_matrix_ij) &
      !$omp private(collocation_matrix_kl, outer_product) &
      !$omp shared(Vijkl, collocation_matrix, Vmunu)
      do i = 1, Norb
        do j = 1, Norb
          collocation_matrix_ij = collocation_matrix(:,(i-1)*Norb+j)
          do k = 1, Norb
            do l = 1, Norb
              collocation_matrix_kl = collocation_matrix(:,(k-1)*Norb+l)
              outer_product = spread(collocation_matrix_ij, 2, nd) * &
                              spread(collocation_matrix_kl, 1, nd)
              Vijkl(i, j, k, l) = sum(Vmunu * outer_product)
            enddo
          enddo
        enddo
      enddo
      !$omp end parallel do

      ! Deallocate arrays
      deallocate(collocation_matrix, collocation_matrix_ij)
      deallocate(collocation_matrix_kl, outer_product)

      endsubroutine computeVijkl
