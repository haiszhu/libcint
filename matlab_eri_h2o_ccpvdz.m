% generate adaptive tree with treefun
%
% 07/28/24 Hai

addpath('./utils/')
addpath('./treefun/')

clear all
order = 7;
% eps = 1e-03; % do not use 1e-03 for dmk call
% eps = 1e-04;
eps = 1e-05;
% eps = 1e-06;
% eps = 1e-07;
func2 = @(x,y,z) cgto2func(x,y,z);
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',[0 0 0; 0 -0.757 0.757;0 0.587 0.587], ...
              'ifcoeffs',false);
f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 
% f = treefun3(func,[-0.5 0.5 -0.5 0.5 -0.5 0.5],order,opts); 
plot(f,func2);

%%% this needs to be changed consistently
Norb = 24; % 
ratio = 0.5/15; % from boxlen to 1

%%% convert to bdmk needed format...
dpars = zeros(1000,1);
ipars = zeros(100,1);
zpars = zeros(10,1);
% dada
dom = ratio*f.domain(:,1)';
scale = (dom(2)-dom(1))/2;
norder = f.n;
ndim = numel(dom(1:end/2));
nhess = ndim*(ndim+1)/2;
% convert to dmk tree format (should not change the tree if opts.balance == true)
iper = 0;
npbox = norder^ndim;
grid = zeros(ndim,norder^2);
nbmax = 2^ndim*f.id(end); % is 2^ndim enough, since some box needs to be refined more than once...
nlmax = f.height(1); % is this ok?
centers0 = ratio/2*(f.domain(1:2:end,:)+f.domain(2:2:end,:));
centers = zeros(ndim,nbmax); centers(:,1:size(centers0,2)) = centers0;
nlevels = f.height(1);
nboxes = f.id(end);
boxsize = 2*scale./2.^(0:nlmax); 
laddr = zeros(2,nlmax+1); 
for ilev=0:nlmax
  laddr(1,ilev+1) = sum(f.level<ilev)+1;
  laddr(2,ilev+1) = sum(f.level<=ilev);
end
ilevel = zeros(nbmax,1); ilevel(1:nboxes) = f.level;
iparent = zeros(nbmax,1); iparent(1:nboxes) = f.parent;
nchild = zeros(nbmax,1); nchild(1:nboxes) = sum(f.children>0);
ichild = zeros(2^ndim,nbmax); ichild(:,1:nboxes) = f.children; ichild(ichild==0) = -1; % -1 convention in fortran tree
% this does not belong to fortran tree, but needed for treefun plot
coeffs = cell(1,nbmax); coeffs(1:nboxes) = f.coeffs(1:nboxes);
% call computecoll to get vol_tree_fix_lr coll info
nnbors0 = zeros(nboxes,1); nbors0 = zeros(3^ndim,nboxes);
[nnbors0,nbors0] = ...
  computecoll(ndim,nlevels,nboxes,laddr,boxsize,...
                  centers(:,1:nboxes),iparent(1:nboxes),nchild(1:nboxes),ichild(:,1:nboxes),iper,...
                  nnbors0,nbors0);
nnbors = zeros(nbmax,1); nbors = zeros(3^ndim,nbmax);
nnbors(1:nboxes) = nnbors0; nbors(:,1:nboxes) = nbors0;
% now fix lr (with coeffs, which is different from)
ndtmp = 1;
[centers,nboxes,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,coeffs] = ...
      vol_tree_fix_lr_treefun(ndim,iper,ndtmp,dpars,zpars,ipars,...
                               norder,npbox,grid,centers,nlevels,nboxes,boxsize,...
                               nbmax,nlmax,laddr,ilevel,iparent,nchild,ichild,nnbors,nbors,...
                               coeffs);
centers = centers(:,1:nboxes);

% after fix lr
ipoly=0;
iperiod=0;
itype=0;
zk=60.0d0;
epstree=opts.tol;
mc=2^ndim;
mnbors=3^ndim;
npbox=norder^ndim;
nleafbox=numel(f.leaves);
boxlen=boxsize(1);
iptype=2;
eta=0.0d0;
rintl=zeros(nlevels+1,1);
iptr=zeros(8,1);
iptr(1) = 1;
iptr(2) = 2*(nlevels+1)+1;
iptr(3) = iptr(2) + nboxes;
iptr(4) = iptr(3) + nboxes;
iptr(5) = iptr(4) + nboxes;
iptr(6) = iptr(5) + mc*nboxes;
iptr(7) = iptr(6) + nboxes;
iptr(8) = iptr(7) + mnbors*nboxes;
ltree = (4 + mc + mnbors) * nboxes + 2 * (nlevels + 1);
%
itree = zeros(ltree,1);
itree(iptr(1):(iptr(1)+2*nlevels+1)) = reshape(laddr,[],1);
itree(iptr(2):(iptr(2)+nboxes-1)) = ilevel(1:nboxes);
itree(iptr(3):(iptr(3)+nboxes-1)) = iparent(1:nboxes);
itree(iptr(4):(iptr(4)+nboxes-1)) = nchild(1:nboxes);
itree(iptr(5):(iptr(5)+2^ndim*nboxes-1)) = reshape(ichild(:,1:nboxes),[],1);
itree(iptr(6):(iptr(6)+nboxes-1)) = nnbors(1:nboxes);
itree(iptr(7):(iptr(7)+3^ndim*nboxes-1)) = reshape(nbors(:,1:nboxes),[],1);
%
src = zeros(ndim,npbox,nboxes);
wts = zeros(npbox,nboxes);
srcleaf = zeros(ndim,npbox,nleafbox);
wtsleaf = zeros(npbox,nleafbox);
xq = zeros(norder, 1);
wtsq = zeros(norder, 1);
umat = zeros(norder, norder);
vmat = zeros(norder, norder);
grid = zeros(ndim, npbox);
wtsb = zeros(npbox, 1);
if ipoly == 0
  [xq, umat, vmat, wtsq] = legeexps(itype, norder, xq, umat, vmat, wtsq);
elseif ipoly == 1
  [xq, umat, vmat, wtsq] = chebexps(itype, norder, xq, umat, vmat, wtsq);
end
xq = xq / 2;
wtsq = wtsq / 2;
grid = mesh3d(xq, norder, xq, norder, xq, norder, grid);
ind = 0;
for iz = 1:norder
  for iy = 1:norder
    for ix = 1:norder
      ind = ind + 1;
      wtsb(ind) = wtsq(ix) * wtsq(iy) * wtsq(iz);
    end
  end
end
% Main loop to compute src on the grid points
jbox = 0;
xyz = zeros(ndim,1);
for ilev = 0:nlevels
  bs = boxsize(ilev+1); % Adjusted for MATLAB indexing (1-based)
  ifirstbox = itree(2 * ilev + 1);
  ilastbox = itree(2 * ilev + 2);
  nbloc = ilastbox - ifirstbox + 1;
  for i = 1:nbloc
    ibox = ifirstbox + i - 1;
    for ell = 1:npbox
      for j = 1:ndim
        xyz(j) = centers(j, ibox) + grid(j, ell) * bs;
        src(j, ell, ibox) = xyz(j);
      end
      wts(ell, ibox) = wtsb(ell) * bs^ndim;
    end
    if itree(iptr(4) + ibox - 1) == 0
      jbox = jbox + 1;
      for ell = 1:npbox
        for j = 1:ndim
          srcleaf(j, ell, jbox) = src(j, ell, ibox);
        end
        wtsleaf(ell, jbox) = wts(ell, ibox);
      end
    end
  end
end

%%% evaluate cGTO on src0 (src are within [-1/2 1/2]^3)
src0 = src/ratio;
func = @(x,y,z) cgtofunc(x,y,z);
fvals0 = squeeze(func(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
fvals0 = permute(fvals0,[3 1 2]);
fvals = fvals0;


%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
nd = Norb*(Norb+1)/2;
phi_kl = zeros(nd,npbox,nboxes); 
idx = 0;
for k = 1:Norb
  for ell = k:Norb
    idx = idx + 1;
    phi_kl(idx,:,:) = fvals(k,:,:).*fvals(ell,:,:);
  end
end
pot = zeros(nd,npbox,nboxes);
nd0 = 10;
nchnk = nd/nd0;
% nchnk = 1;
for jchnk = 1:nchnk
  jidxv = (jchnk-1)*nd0 + (1:nd0);
  phi_kl0 = phi_kl(jidxv,:,:);
  pot0=zeros(nd0,npbox,nboxes);
  grad0=zeros(nd0,ndim,npbox,nboxes);
  hess0=zeros(nd0,nhess,npbox,nboxes);
  tic
  ikernel=1;
  beta=6.0d0;
  eps=opts.tol;
  ifpgh=1;
  ifpghtarg=0;
  ntarg = 100;
  targs=zeros(ndim,ntarg); pote=zeros(nd0,ntarg);
  grade=zeros(nd0,ndim,ntarg); hesse=zeros(nd0,nhess,ntarg);
  timeinfo = zeros(20,1);
  [pot0,grad0,hess0,pote,grade,hesse] = ...
    bdmk_mex(nd0,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
           nboxes,nlevels,ltree,itree,iptr,centers,boxsize,phi_kl0, ...
           ifpgh,pot0,grad0,hess0,ntarg,targs, ...
           ifpghtarg,pote,grade,hesse,timeinfo);
  pot(jidxv,:,:) = pot0;
  time = toc;
end
ntot=npbox*nleafbox+ntarg*ifpghtarg;
X = [' ntotal= ',num2str(ntot)];
disp(X)
pps = ntot/time;
X = [' speed in pps= ', num2str(pps)];
disp(X)


xleaf = srcleaf(1,:,:); xleaf = xleaf(:);
yleaf = srcleaf(2,:,:); yleaf = yleaf(:);
zleaf = srcleaf(3,:,:); zleaf = zleaf(:);
phi_ij_leaf = zeros(nd,npbox,nleafbox); % leaf box information
potleaf = zeros(nd,npbox,nleafbox);
potnonleaf = zeros(nd,npbox,nboxes-nleafbox);
jbox = 0;
jbox2 = 0;
for ilev = 0:nlevels
  bs = boxsize(ilev+1); % Adjusted for MATLAB indexing (1-based)
  ifirstbox = itree(2 * ilev + 1);
  ilastbox = itree(2 * ilev + 2);
  nbloc = ilastbox - ifirstbox + 1;
  for i = 1:nbloc
    ibox = ifirstbox + i - 1;
    if itree(iptr(4) + ibox - 1) == 0
      jbox = jbox + 1;
      for ell = 1:npbox
        for j = 1:nd
          phi_ij_leaf(j, ell, jbox) = phi_kl(j, ell, ibox);
          potleaf(j, ell, jbox) = pot(j, ell, ibox);
        end
      end
    else
      jbox2 = jbox2 + 1;
      for ell = 1:npbox
        for j = 1:nd
          potnonleaf(j, ell, jbox2) = pot(j, ell, ibox);
        end
      end
    end
  end
end

%% compute smooth volume integral
Vijkl0 = zeros(nd,nd);
% ell = 1;
for ell = 1:nd
  % (2*L)^2 * potleaf = (2*L)^2 * \int_{-1/2}^{1/2} (phi_k*phi_l)/|r' - r| dV
  potleaf_ell = potleaf(ell,:,:)/ratio^2; 
  for j = 1:nd
    % (2*L)^3 * \int_{-1/2}^{1/2} phi_{ij} * pot_{kl} dV
    Vijkl0(j,ell) = sum(squeeze(phi_ij_leaf(j,:,:).*potleaf_ell).*wtsleaf(:,:),'all')/ratio^3;
  end
end
% on [-15 15]
% Vijkl1 =  4.738267922428247e+00
% Vijkl2 = -1.600362946654103e+00
% Vijkl3 =  5.961116488992190e-01


keyboard

figure(1),clf,
scatter3(srcleaf(1,1:end/100)',srcleaf(2,1:end/100)',srcleaf(3,1:end/100)',2,potleaf_ell(1:end/100))

% turn back into 4-tensor and save
triuidx = 1:Norb^2;
triflag = true(Norb,Norb);
triflag = triu(triflag);
triuidx = triuidx(triflag(:));
Vijkltmp = zeros(Norb,Norb,nd);
for kl = 1:nd
  tmpidx = 0;
  for i = 1:Norb
    for j = i:Norb
      tmpidx = tmpidx + 1;
      Vijkltmp(i,j,kl) = Vijkl0(tmpidx,kl);
    end
  end
  for i = 1:Norb
    for j = 1:i-1
      Vijkltmp(i,j,kl) = Vijkltmp(j,i,kl);
    end
  end
end
Vijkl = zeros(Norb,Norb,Norb,Norb);
tmpidx = 0;
for k = 1:Norb
  for l = k:Norb
    tmpidx = tmpidx + 1;
    Vijkl(:,:,k,l) = Vijkltmp(:,:,tmpidx);
  end
end
for k = 1:Norb
  for l = 1:k-1
    Vijkl(:,:,k,l) = Vijkl(:,:,l,k);
  end
end
save('ERI_h2o_ccpvdz.mat','Vijkl')    
h5create("ERI_h2o_ccpvdz.h5","/DS1",[Norb Norb Norb Norb])
h5write("ERI_h2o_ccpvdz.h5","/DS1",Vijkl)
h5disp("ERI_h2o_ccpvdz.h5")