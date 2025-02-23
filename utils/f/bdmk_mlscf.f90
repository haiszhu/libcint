      ! export OPENBLAS_NUM_THREADS=1
      ! export OMP_NUM_THREADS=1 
      ! gfortran -o int2-bdmk-mlscf -I/opt/homebrew/include -L/opt/homebrew/lib -lhdf5_fortran -lhdf5 -fopenmp bdmk_mlscf.f90 bdmk_wrap.o bdmk.o boxfgt_md.o ./src/besseljs3d.o ./src/hank103.o ./src/legeexps.o ./src/chebexps.o ./src/prini_new.o ./src/fmmcommon2d.o ./src/lapack_f77.o ./src/cumsum.o ./src/hkrand.o ./src/dlaran.o ./src/voltab2d.o ./src/voltab3d.o ./src/polytens.o ./src/tree_data_routs.o ./src/tensor_prod_routs.o ./src/pts_tree.o ./src/tree_routs.o ./src/tree_vol_coeffs.o ./src/dmk_routs.o ./src/get_sognodes.o ./src/l2dsognodes.o ./src/l3dsognodes.o ./src/sl3dsognodes.o ./src/y2dsognodes.o ./src/y3dsognodes.o ./src/bdmk_local_tables.o ./src/bdmk_local.o ./src/bdmk_pwterms.o ./src/bdmk_pwrouts.o -L. /Users/hzhu/git/OpenBLAS/libopenblas_vortex-r0.3.29.a -lpthread -lm 
      ! gfortran -O3 -o int2-bdmk-mlscf -fopenmp bdmk_mlscf.f90 bdmk_wrap.o bdmk.o boxfgt_md.o ./src/besseljs3d.o ./src/hank103.o ./src/legeexps.o ./src/chebexps.o ./src/prini_new.o ./src/fmmcommon2d.o ./src/lapack_f77.o ./src/cumsum.o ./src/hkrand.o ./src/dlaran.o ./src/voltab2d.o ./src/voltab3d.o ./src/polytens.o ./src/tree_data_routs.o ./src/tensor_prod_routs.o ./src/pts_tree.o ./src/tree_routs.o ./src/tree_vol_coeffs.o ./src/dmk_routs.o ./src/get_sognodes.o ./src/l2dsognodes.o ./src/l3dsognodes.o ./src/sl3dsognodes.o ./src/y2dsognodes.o ./src/y3dsognodes.o ./src/bdmk_local_tables.o ./src/bdmk_local.o ./src/bdmk_pwterms.o ./src/bdmk_pwrouts.o -L. /home/hai/git/OpenBLAS/libopenblas.a -lpthread -lm 
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
      real *8 eps, beta
      integer nd, ndim, i, j, k, ios
      integer ikernel, nboxes, nlevels, norder, npbox, ipoly, ltree
      integer iptr(8)      
      integer, allocatable :: itree(:)
      real *8, allocatable :: centers(:,:)
      real *8, allocatable :: boxsize(:)
      real *8, allocatable :: fvals(:,:,:)
      real *8, allocatable :: pot(:,:,:)
      ! 
      real *8 ndr, ndimr
      real *8 ikernelr,nboxesr,nlevelsr,norderr,npboxr,ipolyr,ltreer
      real *8 iptrr(8)
      real *8, allocatable :: itreer(:)
      real *8, allocatable :: phi_kl(:,:)
      real *8 h5diff
      ! 
      character(len=100) :: filename1, filename2
      character(len=100) :: treefun_h5_fn
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
      timeinfo = 0.0d0
      ! 
            
      call h5open_f(error) ! Initialize HDF5
      ! read treefun data
      call h5fopen_f("treefun_h2o_cc-pvdz.h5", H5F_ACC_RDONLY_F, & 
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
      ! treefun_file_id_rd, treefun_dataset_id_rd

      call h5dclose_f(treefun_dataset_id_rd, error) ! close original dataset and file

      ! read isdf data
      call h5fopen_f("isdf_1e-3.h5", H5F_ACC_RDONLY_F, & 
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
      ! call cpu_time(t1)
      t1 = omp_get_wtime()
      call bdmk_wrap(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, &
          nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals, &
          pot)
      ! call cpu_time(t2)
      t2 = omp_get_wtime()
      timeinfo(1) = timeinfo(1) + (t2-t1)

      PRINT *, 'bdmk time is', timeinfo(1), ' seconds.'

      ! write
      dims_wt = (/nd, npbox, nboxes/)
      call h5fcreate_f("bdmk_1e-3.h5", H5F_ACC_TRUNC_F,isdf_file_id_wt,error) ! Create a new HDF5 file
      call h5screate_simple_f(3, dims_wt, isdf_dataspace_id_wt, error) ! Create a dataspace for the new dataset
      call h5dcreate_f(isdf_file_id_wt, "pot", H5T_NATIVE_DOUBLE, &
                    isdf_dataspace_id_wt, isdf_dataset_id_wt, error)
      call h5dwrite_f(isdf_dataset_id_wt, H5T_NATIVE_DOUBLE, pot, &
                    dims_wt, error) 
      if (error /= 0) then
        print*, "Error occurred while writing HDF5 file."
      endif
      call h5dclose_f(isdf_dataset_id_wt, error) ! Close dataset and file
      call h5sclose_f(isdf_dataspace_id_wt, error)
      call h5fclose_f(isdf_file_id_wt, error)
      call h5close_f(error) 

      deallocate(itree)
      deallocate(centers)
      deallocate(boxsize)
      deallocate(interpolating_vectors)
      deallocate(fvals)
      deallocate(pot)
      deallocate(DS1)

      end program bdmk_mlscf