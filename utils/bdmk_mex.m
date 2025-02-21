function [pot,grad,hess,pote,grade,hesse] = bdmk_mex(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals,ifpgh,pot,grad,hess,ntarg,targs,ifpghtarg,pote,grade,hesse,tottimeinfo)
nlevelsp1 = nlevels+1;
nhess = ndim*(ndim+1)/2;
npboxnboxes = npbox*nboxes;
ndimnpboxnboxes = ndim*npbox*nboxes;
nhessnpboxnboxes = nhess*npboxnboxes;
fvals = reshape(fvals,[nd npboxnboxes]);
pot = reshape(pot,[nd npboxnboxes]);
grad = reshape(grad,[nd ndimnpboxnboxes]);
hess = reshape(hess,[nd nhessnpboxnboxes]);
ndimntarg = ndim*ntarg;
nhessntarg = nhess*ntarg;
grade = reshape(grade,[nd ndimntarg]);
hesse = reshape(hesse,[nd nhessntarg]);
mex_id_ = 'bdmk(i int[x], i int[x], i double[x], i int[x], i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], i int[x], i double[xx], i int[x], io double[xx], io double[xx], io double[xx], io double[x])';
[pot, grad, hess, pote, grade, hesse, tottimeinfo] = bdmk_module(mex_id_, nd, ndim, eps, ikernel, beta, ipoly, norder, npbox, nboxes, nlevels, ltree, itree, iptr, centers, boxsize, fvals, ifpgh, pot, grad, hess, ntarg, targs, ifpghtarg, pote, grade, hesse, tottimeinfo, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ltree, 8, ndim, nboxes, nlevelsp1, nd, npboxnboxes, 1, nd, npboxnboxes, nd, ndimnpboxnboxes, nd, nhessnpboxnboxes, 1, ndim, ntarg, 1, nd, ntarg, nd, ndimntarg, nd, nhessntarg, 20);
pot = reshape(pot,[nd npbox nboxes]);
grad = reshape(grad,[nd ndim npbox nboxes]);
hess = reshape(hess,[nd nhess npbox nboxes]);
grade = reshape(grade,[nd ndim ntarg]);
hesse = reshape(hesse,[nd nhess ntarg]);
end