function pot = bdmk_warp(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                         nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvals, ...
                         pot)
nhess = ndim*(ndim+1)/2;
%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
nd0 = 1;
disp("=========Start BDMK=======");
disp("nchnk is : " + nd );
for jchnk = 1:nd
  fvalsj = fvals(jchnk,:,:);
  pot0=zeros(nd0,npbox,nboxes);
  grad0=zeros(nd0,ndim,npbox,nboxes);
  hess0=zeros(nd0,nhess,npbox,nboxes);
  tic
  ifpgh=1;
  ifpghtarg=0;
  ntarg = 2;
  targs=zeros(ndim,ntarg); pote=zeros(nd0,ntarg);
  grade=zeros(nd0,ndim,ntarg); hesse=zeros(nd0,nhess,ntarg);
  timeinfo = zeros(20,1);
  [pot0,grad0,hess0,pote,grade,hesse] = ...
    bdmk_mex(nd0,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
           nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvalsj, ...
           ifpgh,pot0,grad0,hess0,ntarg,targs, ...
           ifpghtarg,pote,grade,hesse,timeinfo);
  pot(jchnk,:,:) = pot0;
  disp("jchnk is : " + jchnk + " ( total nchnk is : " + nd + " )");
end
disp("=========End BDMK=======");

end