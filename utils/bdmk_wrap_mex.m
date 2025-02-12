function [pot] = bdmk_wrap_mex(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,pot)
nlevelsp1 = nlevels+1;
npboxnboxes = npbox*nboxes;
fvals = reshape(fvals,[nd npboxnboxes]);
mex_id_ = 'bdmk_wrap(i int[x], i int[x], i double[x], i int[x], i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i double[xx], io double[xx])';
[pot] = bdmk_module(mex_id_, nd, ndim, eps, ikernel, beta, ipoly, norder, npbox, nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals, pot, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ltree, 8, ndim, nboxes, nlevelsp1, nd, npboxnboxes, nd, npboxnboxes);
end