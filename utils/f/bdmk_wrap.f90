      subroutine bdmk_wrap(nd,ndim,eps,ikernel,beta,ipoly,norder, &
          npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,  &
          fvals,pot)
      use omp_lib
      implicit none
      real *8 eps, beta
      integer nd, ndim
      integer ikernel, nboxes, nlevels, norder, npbox, ipoly, ltree
      integer itree(ltree), iptr(8)
      real *8 centers(ndim,nboxes)
      real *8 boxsize(0:nlevels)
      real *8 fvals(nd,npbox,nboxes)
      real *8 pot(nd,npbox,nboxes)
      ! 
      integer jchnk, nd0, nhess, ifpgh, ifpghtarg, ntarg
      real *8, allocatable :: fvalsj(:,:,:)
      real *8, allocatable :: pot0(:,:,:),grad0(:,:,:,:),hess0(:,:,:,:)
      real *8, allocatable :: targs(:,:),pote(:,:),grade(:,:,:)
      real *8, allocatable :: hesse(:,:,:),timeinfo(:)
      ! 
      nd0 = 1
      nhess = ndim*(ndim+1)/2
      ifpgh = 1
      ifpghtarg = 0
      ntarg = 2
      ! 
      !$omp parallel default(none) &
      !$omp& shared(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox) & 
      !$omp& shared(nboxes,nlevels,ltree,itree,iptr,centers,boxsize) &
      !$omp& shared(fvals,pot,nd0,nhess,ifpgh,ifpghtarg,ntarg) &
      !$omp& private(jchnk,fvalsj,pot0,grad0,hess0,targs,pote) &
      !$omp& private(grade,hesse,timeinfo)
      ! 
      allocate(fvalsj(nd0, npbox, nboxes))
      allocate(pot0(nd0, npbox, nboxes))
      allocate(grad0(nd0, ndim, npbox, nboxes))
      allocate(hess0(nd0, nhess, npbox, nboxes))
      allocate(targs(ndim, ntarg))
      allocate(pote(nd0, ntarg))
      allocate(grade(nd0, ndim, ntarg))
      allocate(hesse(nd0, nhess, ntarg))
      allocate(timeinfo(20))
      !$omp do 
      do jchnk = 1,nd
        ! print *, 'I am at jchnk = ', jchnk
        ! 
        fvalsj = 0.0d0
        pot0   = 0.0d0
        grad0  = 0.0d0
        hess0  = 0.0d0
        targs  = 0.0d0
        pote   = 0.0d0
        grade  = 0.0d0
        hesse  = 0.0d0
        timeinfo = 0.0d0
        ! 
        fvalsj(1,:,:) = fvals(jchnk, :, :)
        ! 
        call bdmk(nd0,ndim,eps,ikernel,beta,ipoly,norder,npbox, &
          nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvalsj, &
          ifpgh,pot0,grad0,hess0,ntarg,targs, &
          ifpghtarg,pote,grade,hesse,timeinfo)
        ! 
        pot(jchnk, :, :) = pot0(1, :, :)
        ! print *, 'I am done with jchnk = ', jchnk
      enddo
      !$omp end do
      ! 
      deallocate(fvalsj,pot0,grad0,hess0)
      deallocate(targs,pote,grade,hesse,timeinfo)
      !$omp end parallel
      
      end subroutine bdmk_wrap