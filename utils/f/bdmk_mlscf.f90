      ! export OPENBLAS_NUM_THREADS=1
      ! export OMP_NUM_THREADS=1 
      ! on max:   gfortran -O3 -o int2-bdmk-mlscf -I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5_fortran -lhdf5 -fopenmp bdmk_mlscf.f90 bdmk_wrap.o bdmk.o boxfgt_md.o ./src/besseljs3d.o ./src/hank103.o ./src/legeexps.o ./src/chebexps.o ./src/prini_new.o ./src/fmmcommon2d.o ./src/lapack_f77.o ./src/cumsum.o ./src/hkrand.o ./src/dlaran.o ./src/voltab2d.o ./src/voltab3d.o ./src/polytens.o ./src/tree_data_routs.o ./src/tensor_prod_routs.o ./src/pts_tree.o ./src/tree_routs.o ./src/tree_vol_coeffs.o ./src/dmk_routs.o ./src/get_sognodes.o ./src/l2dsognodes.o ./src/l3dsognodes.o ./src/sl3dsognodes.o ./src/y2dsognodes.o ./src/y3dsognodes.o ./src/bdmk_local_tables.o ./src/bdmk_local.o ./src/bdmk_pwterms.o ./src/bdmk_pwrouts.o -L. /Users/hzhu/git/OpenBLAS/libopenblas_vortex-r0.3.29.a -lpthread -lm 
      ! on linux: gfortran -O3 -o int2-bdmk-mlscf -fopenmp bdmk_mlscf.f90 bdmk_wrap.o bdmk.o boxfgt_md.o ./src/besseljs3d.o ./src/hank103.o ./src/legeexps.o ./src/chebexps.o ./src/prini_new.o ./src/fmmcommon2d.o ./src/lapack_f77.o ./src/cumsum.o ./src/hkrand.o ./src/dlaran.o ./src/voltab2d.o ./src/voltab3d.o ./src/polytens.o ./src/tree_data_routs.o ./src/tensor_prod_routs.o ./src/pts_tree.o ./src/tree_routs.o ./src/tree_vol_coeffs.o ./src/dmk_routs.o ./src/get_sognodes.o ./src/l2dsognodes.o ./src/l3dsognodes.o ./src/sl3dsognodes.o ./src/y2dsognodes.o ./src/y3dsognodes.o ./src/bdmk_local_tables.o ./src/bdmk_local.o ./src/bdmk_pwterms.o ./src/bdmk_pwrouts.o -L. /home/hai/git/OpenBLAS/libopenblas.a -lpthread -lm -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 
      ! 
      ! on linux: gfortran -O3 -o int2-bdmk-mlscf -fopenmp bdmk_mlscf.f90 bdmk_wrap.o bdmk.o boxfgt_md.o besseljs3d.o hank103.o legeexps.o chebexps.o prini_new.o fmmcommon2d.o lapack_f77.o cumsum.o hkrand.o dlaran.o voltab2d.o voltab3d.o polytens.o tree_data_routs.o tensor_prod_routs.o pts_tree.o tree_routs.o tree_vol_coeffs.o dmk_routs.o get_sognodes.o l2dsognodes.o l3dsognodes.o sl3dsognodes.o y2dsognodes.o y3dsognodes.o bdmk_local_tables.o bdmk_local.o bdmk_pwterms.o bdmk_pwrouts.o -L. /home/hai/git/OpenBLAS/libopenblas.a -lpthread -lm -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 
      program bdmk_mlscf
      ! 
      ! 
      use omp_lib 
      use hdf5
      implicit none
      ! isdf .h5 variables
      ! treefun
      integer(HID_T) :: treefun_file_id_rd, treefun_dataset_id_rd
      integer(HID_T) :: treefun_dataspace_id
      integer(HSIZE_T), dimension(1) :: scalar_dims
      integer(HSIZE_T), dimension(2) :: vector_dims
      integer(HSIZE_T), dimension(3) :: tensor_dims
      ! integer(HSIZE_T), dimension(1) :: dims
      ! isdf
      integer(HID_T) :: isdf_file_id_rd, isdf_dataset_id_rd
      integer(HID_T) :: dataspace_id
      integer(HID_T) :: isdf_file_id_wt, isdf_dataset_id_wt 
      integer(HID_T) :: isdf_dataspace_id_wt
      integer(HSIZE_T), dimension(3) :: dims_rd, dims_wt
      integer :: error
      real *8, allocatable :: interpolating_vectors(:,:,:)
      integer *8 :: maxdims(3) ! this is a bit awkward...
      ! bdmk variables
      real *8 eps, beta, ratio
      integer nd, ndim, i, j, k, ios
      integer ikernel, nboxes, nlevels, norder, npbox, ipoly, ltree
      integer iptr(8)      
      integer, allocatable :: itree(:)
      real *8, allocatable :: centers(:,:)
      real *8, allocatable :: boxsize(:)
      real *8, allocatable :: fvals(:,:,:)
      real *8, allocatable :: pot(:,:,:)
      integer nleafbox
      real *8, allocatable :: wtsleaf(:,:)
      ! 
      real *8 ndr, ndimr
      real *8 ikernelr,nboxesr,nlevelsr,norderr,npboxr,ipolyr,ltreer
      real *8 iptrr(8)
      real *8, allocatable :: itreer(:)
      real *8, allocatable :: phi_kl(:,:)
      real *8 h5diff
      ! 
      character(len=100) :: treefun_filename
      character(len=100) :: isdf_filename
      character(len=100) :: output_filename
      integer :: num_args
      ! 
      integer ndims
      ! integer *8 :: dims(1)
      ! integer(HSIZE_T), dimension(1) :: dims
      ! integer(8), dimension(1) :: dims
      ! integer(4), dimension(1) :: dims
      integer(kind=HSIZE_T), dimension(1) :: dims
      real *8, allocatable :: DS1(:,:,:)
      ! 
      real *8 :: timeinfo(20)
      real *8 :: t1, t2
      ! 
      real *8 :: bs
      real *8, allocatable :: potleaf(:,:,:), fvalsleaf(:,:,:)
      integer :: jbox, ibox, ilev, ifirstbox, ilastbox, nbloc, ell
      ! Vmunu
      real *8, allocatable :: potleafrs(:,:), fvalsleafrs(:,:)
      real *8, allocatable :: wtsleafrs(:), Vmunu2(:,:)
      real *8, allocatable :: Vmunu(:,:)
      real *8 :: alpha
      integer :: m, n, lda, ldb, ldc
      ! 
      timeinfo = 0.0d0
      !
      num_args = command_argument_count()
      !
      if (num_args >= 3) then
        call get_command_argument(1, treefun_filename)
        call get_command_argument(2, isdf_filename)
        call get_command_argument(3, output_filename)
        print *, 'Using the following filenames:'
        print *, 'treefun_filename: ', trim(treefun_filename)
        print *, 'isdf_filename:    ', trim(isdf_filename)
        print *, 'output_filename:  ', trim(output_filename)
      else
        treefun_filename = "treefun_h2o_cc-pvdz.h5"
        isdf_filename = "isdf_1e-3.h5"
        output_filename = "bdmk_1e-3.h5"
      endif
      !      
      call h5open_f(error) ! Initialize HDF5
      ! read treefun data
      call h5fopen_f(treefun_filename, H5F_ACC_RDONLY_F, & 
                    treefun_file_id_rd,error) ! Open the HDF5 file
      if (error /= 0) then
        print*, "Error opening HDF5 file: "
        stop
      endif
      scalar_dims = (/1/)
      ! read ndim
      call h5dopen_f(treefun_file_id_rd, "/ndim", &
                    treefun_dataset_id_rd, error) ! Open the dataset
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      ndim, scalar_dims, error) ! Read the data into the Fortran array
      ! read treefun_eps
      call h5dopen_f(treefun_file_id_rd, "/eps", &
                    treefun_dataset_id_rd, error) 
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, &
                      eps, scalar_dims, error) 
      ! read ikernel
      call h5dopen_f(treefun_file_id_rd, "/ikernel", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      ikernel, scalar_dims, error)
      ! read beta
      call h5dopen_f(treefun_file_id_rd, "/beta", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, &
                      beta, scalar_dims, error)
      ! read ipoly
      call h5dopen_f(treefun_file_id_rd, "/ipoly", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      ipoly, scalar_dims, error)
      ! read norder
      call h5dopen_f(treefun_file_id_rd, "/norder", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      norder, scalar_dims, error)
      ! read npbox
      call h5dopen_f(treefun_file_id_rd, "/npbox", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      npbox, scalar_dims, error)
      ! read nboxes
      call h5dopen_f(treefun_file_id_rd, "/nboxes", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      nboxes, scalar_dims, error)
      ! read nlevels
      call h5dopen_f(treefun_file_id_rd, "/nlevels", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      nlevels, scalar_dims, error)
      ! read ltree
      call h5dopen_f(treefun_file_id_rd, "/ltree", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      ltree, scalar_dims, error)
      ! read itree
      allocate(itree(ltree))
      scalar_dims = (/ltree/)
      call h5dopen_f(treefun_file_id_rd, "/itree", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      itree, scalar_dims, error)
      ! read iptr
      scalar_dims = (/8/)
      call h5dopen_f(treefun_file_id_rd, "/iptr", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      iptr, scalar_dims, error)
      ! read centers(:)
      allocate(centers(ndim,nboxes))
      vector_dims = (/ndim,nboxes/)
      call h5dopen_f(treefun_file_id_rd, "/centers", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, & 
                      centers, vector_dims, error)
      ! read boxsize
      allocate(boxsize(0:nlevels))
      scalar_dims = (/nlevels+1/)
      call h5dopen_f(treefun_file_id_rd, "/boxsize", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, &
                      boxsize, scalar_dims, error)
      ! read DS1, i.e. src
      allocate(DS1(ndim,npbox,nboxes))
      tensor_dims = (/ndim,npbox,nboxes/)
      call h5dopen_f(treefun_file_id_rd, "/DS1", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, &
                      DS1, tensor_dims, error)
      ! read nleafbox
      scalar_dims = (/1/)
      call h5dopen_f(treefun_file_id_rd, "/nleafbox", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_INTEGER, &
                      nleafbox, scalar_dims, error)
      ! read wtsleaf
      allocate(wtsleaf(npbox,nleafbox))
      vector_dims = (/npbox,nleafbox/)
      call h5dopen_f(treefun_file_id_rd, "/wtsleaf", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, & 
                      wtsleaf, vector_dims, error)
      ! read ratio
      call h5dopen_f(treefun_file_id_rd, "/ratio", &
                    treefun_dataset_id_rd, error)
      call h5dread_f(treefun_dataset_id_rd, H5T_NATIVE_DOUBLE, & 
                      ratio, scalar_dims, error)

      print *, "Size of HSIZE_T:", kind(HSIZE_T)
      print *, 'ndim should be: ', ndim
      print *, 'ndim read from .h5 is: ', ndim
      print *, 'ndim is ', ndim
      print *, 'eps is ', eps
      print *, 'ikernel is ', ikernel
      print *, 'beta is ', beta
      print *, 'ipoly is ', ipoly
      print *, 'norder is ', norder
      print *, 'npbox is ', npbox
      print *, 'nboxes is ', nboxes
      print *, 'ltree is ', ltree
      print *, 'itree is ', itree(1:10)
      print *, 'iptr is ', iptr
      print *, 'center is ', centers(1,1), centers(2,1), centers(3,1), &
                             centers(1,2), centers(2,2), centers(3,2), &
                             centers(1,3), centers(2,3), centers(3,3)
      print *, 'boxsize is ', boxsize(:)
      print *, 'nleafbox is ', nleafbox
      print *, 'wtsleaf is ', wtsleaf(1:10,1)
      print *, 'ratio is ', ratio
      ! treefun_file_id_rd, treefun_dataset_id_rd

      call h5dclose_f(treefun_dataset_id_rd, error) ! close original dataset and file

      ! read isdf data
      call h5fopen_f(isdf_filename, H5F_ACC_RDONLY_F, & 
                    isdf_file_id_rd,error) ! Open the HDF5 file
      ! call h5fopen_f("isdf_1e-3.h5", H5F_ACC_RDWR_F, isdf_file_id_rd, error) ! Open the HDF5 file
      call h5dopen_f(isdf_file_id_rd, "/interpolating_vectors", &
                    isdf_dataset_id_rd, error) ! Open the dataset
      call h5dget_space_f(isdf_dataset_id_rd, dataspace_id, error)
      call h5sget_simple_extent_dims_f(dataspace_id, dims_rd, maxdims, error)
      ! print *, "Dimensions:", dims_rd
      allocate(interpolating_vectors(dims_rd(1),dims_rd(2),dims_rd(3))) ! Allocate array to store data
      call h5dread_f(isdf_dataset_id_rd, H5T_NATIVE_DOUBLE, &
                      interpolating_vectors, dims_rd, error) ! Read the data into the Fortran array
      call h5dclose_f(isdf_dataset_id_rd, error) ! close original dataset and file
      call h5fclose_f(isdf_file_id_rd, error)

      ! bdmk input format
      nd = dims_rd(3)
      allocate(fvals(nd,npbox,nboxes))
      do k=1,nd
        do j=1,nboxes
          do i=1,npbox
            fvals(k,i,j) = interpolating_vectors(1,i+(j-1)*npbox,k)
          enddo
        enddo
      enddo

      allocate(pot(nd,npbox,nboxes))
      PRINT *, '=========Start BDMK======='
      ! call cpu_time(t1)
      t1 = omp_get_wtime()
      call bdmk_wrap(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, &
          nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals, &
          pot)
      ! call cpu_time(t2)
      t2 = omp_get_wtime()
      timeinfo(1) = (t2-t1)
      PRINT *, '=========End BDMK======='
      PRINT *, 'bdmk time is', timeinfo(1), ' seconds.'

      ! 
      PRINT *, '=========Start potleaf======='
      allocate(potleaf(nd, npbox, nleafbox))
      allocate(fvalsleaf(nd, npbox, nleafbox))
      jbox = 0
      t1 = omp_get_wtime()
      do ilev = 0, nlevels
        bs = boxsize(ilev)
        ifirstbox = itree(2*ilev+1)
        ilastbox = itree(2*ilev+2)
        nbloc = ilastbox-ifirstbox+1
        do i = 1,nbloc
          ibox = ifirstbox+i-1
          if (itree(iptr(4)+ibox-1) == 0) then
            jbox = jbox + 1
            do ell = 1,npbox
              do j = 1,nd
                fvalsleaf(j, ell, jbox) = fvals(j, ell, ibox)
                potleaf(j, ell, jbox) = pot(j, ell, ibox)
              enddo
            enddo
          endif
        enddo
      enddo
      t2 = omp_get_wtime()
      timeinfo(2) = (t2-t1)
      PRINT *, '=========End potleaf======='
      PRINT *, 'potleaf time is', timeinfo(2), ' seconds.'

      ! ! 
      ! PRINT *, '=========Start Vmunu======='
      ! allocate(Vmunu(nd, nd))
      ! Vmunu = 0.0d0
      ! t1 = omp_get_wtime()
      ! do ell = 1,nd
      !   do j = 1,nd
      !     Vmunu(j,ell) = 0.0d0
      !     do i = 1,npbox
      !       do k = 1,nleafbox
      !         Vmunu(j, ell) = Vmunu(j, ell) + fvalsleaf(j, i, k) * &
      !             (potleaf(ell, i, k)/ratio**2) * wtsleaf(i, k)
      !       enddo
      !     enddo
      !   enddo
      ! enddo
      ! Vmunu = Vmunu/ratio**3
      ! t2 = omp_get_wtime()
      ! timeinfo(3) = (t2-t1)
      ! PRINT *, '=========End Vmunu======='
      ! PRINT *, 'Vmunu time is ', timeinfo(3), ' seconds.'

      PRINT *, '=========Start Vmunu (matmul)======='
      t1 = omp_get_wtime()
      allocate(Vmunu2(nd, nd))
      allocate(potleafrs(nd, npbox * nleafbox))
      allocate(fvalsleafrs(nd, npbox * nleafbox))
      allocate(wtsleafrs(npbox * nleafbox))
      potleafrs = reshape(potleaf, [nd, npbox * nleafbox])*(1.0d0/ratio**2)
      fvalsleafrs = reshape(fvalsleaf, [nd, npbox * nleafbox])
      wtsleafrs = reshape(wtsleaf, [npbox * nleafbox])
      fvalsleafrs = fvalsleafrs * spread(wtsleafrs, 1, nd)
      ! ! Vmunu2 = matmul(fvalsleafrs, transpose(potleafrs))
      ! !$omp parallel default(none) &
      ! !$omp& shared(nd, potleafrs, fvalsleafrs, Vmunu2) &
      ! !$omp& private(ell) 
      ! !$omp do
      ! do ell = 1,nd
      !   Vmunu2(:, ell) = matmul(fvalsleafrs, potleafrs(ell, :))
      ! enddo
      ! !$omp end do
      ! !$omp end parallel

      alpha = 1.0d0
      beta = 0.0d0
      m = nd
      n = 1  ! Single column output
      k = npbox * nleafbox  ! Inner dimension
      lda = nd
      ldb = k
      ldc = nd
      !$omp parallel do default(none) &
      !$omp& shared(nd, npbox, nleafbox, fvalsleafrs, potleafrs, Vmunu2, alpha, beta, m, n, k, lda, ldb, ldc)
      do ell = 1, nd
          call dgemm('N', 'N', m, n, k, alpha, fvalsleafrs, lda, potleafrs(ell, :), ldb, beta, Vmunu2(:, ell), ldc)
      end do
      !$omp end parallel do
      
      Vmunu2 =  (1.0d0/ratio**3)*Vmunu2
      Vmunu = Vmunu2
      t2 = omp_get_wtime()
      timeinfo(3) = (t2-t1)
      PRINT *, '=========End Vmunu======='
      PRINT *, 'Vmunu time is ', timeinfo(3), ' seconds.'
      ! PRINT *, 'max diff in Vmunu is ', maxval(Vmunu - Vmunu2)

      ! write
      call h5fcreate_f(output_filename, H5F_ACC_TRUNC_F,isdf_file_id_wt,error) ! Create a new HDF5 file
      ! ! write pot
      ! dims_wt = (/nd, npbox, nboxes/)
      ! call h5screate_simple_f(3, dims_wt, isdf_dataspace_id_wt, error) ! Create a dataspace for the new dataset
      ! call h5dcreate_f(isdf_file_id_wt, "pot", H5T_NATIVE_DOUBLE, &
      !               isdf_dataspace_id_wt, isdf_dataset_id_wt, error)
      ! call h5dwrite_f(isdf_dataset_id_wt, H5T_NATIVE_DOUBLE, pot, &
      !               dims_wt, error) 
      ! if (error == 0) then
      !   print*, "Write variable pot to HDF5 file."
      ! endif
      ! ! write potleaf
      ! dims_wt = (/nd, npbox, nleafbox/)
      ! call h5screate_simple_f(3, dims_wt, isdf_dataspace_id_wt, error) ! Create a dataspace for the new dataset
      ! call h5dcreate_f(isdf_file_id_wt, "potleaf", H5T_NATIVE_DOUBLE, &
      !               isdf_dataspace_id_wt, isdf_dataset_id_wt, error)
      ! call h5dwrite_f(isdf_dataset_id_wt, H5T_NATIVE_DOUBLE, potleaf, &
      !               dims_wt, error) 
      ! if (error == 0) then
      !   print*, "Write variable potleaf to HDF5 file."
      ! endif
      ! write Vmunu
      vector_dims = (/nd,nd/)
      call h5screate_simple_f(2,vector_dims,isdf_dataspace_id_wt,error)
      call h5dcreate_f(isdf_file_id_wt, "Vmunu", H5T_NATIVE_DOUBLE, &
                    isdf_dataspace_id_wt, isdf_dataset_id_wt, error)
      call h5dwrite_f(isdf_dataset_id_wt, H5T_NATIVE_DOUBLE, Vmunu, &
                    dims_wt, error)
      if (error == 0) then
        print*, "Write variable Vmunu to HDF5 file."
      endif

      ! if (error /= 0) then
      !   print*, "Error occurred while writing HDF5 file."
      ! endif
      call h5dclose_f(isdf_dataset_id_wt, error) ! Close dataset and file
      call h5sclose_f(isdf_dataspace_id_wt, error)
      call h5fclose_f(isdf_file_id_wt, error)
      call h5close_f(error) 

      deallocate(itree)
      deallocate(centers)
      deallocate(boxsize)
      deallocate(interpolating_vectors)
      deallocate(fvals)
      deallocate(pot,potleaf,fvalsleaf)
      deallocate(DS1)
      deallocate(wtsleaf)
      deallocate(Vmunu)
      deallocate(Vmunu2)
      deallocate(potleafrs)
      deallocate(fvalsleafrs)
      deallocate(wtsleafrs)

      end program bdmk_mlscf