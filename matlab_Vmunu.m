% dummy test 2
% 10/21/24 Hai

addpath('./utils/')
addpath('./treefun/')
addpath('./data')

clear all
order = 8;
eps = 1e-06; 

%%% resolve tree on cgto^2
func2 = @(x,y,z) cgto2func(x,y,z);
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',[0 0 0; 0 -0.757 0.757;0 0.587 0.587], ...
              'ifcoeffs',false);
f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 
% f = treefun3(func,[-0.5 0.5 -0.5 0.5 -0.5 0.5],order,opts); 
plot(f,func2);

%%% treefun to bdmk
Norb = 24; % 
ndim = 3;
ratio = 0.5/15; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

r = h5read('src_h2o_ccpvdz.h5', '/DS1');
npts = numel(r(1,:));
r2 = reshape(r,3,[]);
diff = abs(reshape(src/ratio,3,[]) - r2);

%%% eval cgto
src0 = src/ratio;
func = @(x,y,z) cgtofunc(x,y,z);
fvals0 = squeeze(func(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
fvals0 = permute(fvals0,[3 1 2]);
fvals = fvals0;

%%% compute V_ijkl
nd = Norb*(Norb+1)/2;
ikernel = 1;
beta = 6.0d0;
% 
% Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
%                 ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
%                 nboxes,nlevels,ltree,itree,iptr,centers,boxsize);

%
idfname = 'isdf_1e-9.h5'; % 0.0026
% idfname = 'isdf_1e-4.h5'; % 0.4476
info = h5info(idfname);
Np = h5read(idfname, '/Np');
collocation_matrix = h5read(idfname, '/collocation_matrix');
interpolating_points = h5read(idfname, '/interpolating_points');
interpolating_vectors = h5read(idfname, '/interpolating_vectors');
kpts = h5read(idfname, '/kpts');
nkpts_ibz = h5read(idfname, '/nkpts_ibz');
nqpts_ibz = h5read(idfname, '/nqpts_ibz');
qpts = h5read(idfname, '/qpts');
  
%
collocation_matrix = squeeze(collocation_matrix(1,:,:));
interpolating_vectors = squeeze(interpolating_vectors(1,:,:));
rho_ij_id = zeros(npts,Norb^2);
for i=1:Norb
  for j=1:Norb
    rho_ij_id(:,(i-1)*Norb+j) = interpolating_vectors...
                                  *(collocation_matrix(i,:).*collocation_matrix(j,:))';

  end
end

func = @(x,y,z) cgtofunc(x,y,z);
fvals0 = func(r2(1,:),r2(2,:),r2(3,:));
fvals0 = reshape(fvals0,[npts Norb]);
fvals_ij = zeros(npts,Norb^2);
for i=1:Norb
  for j=1:Norb
    fvals_ij(:,(i-1)*Norb+j) = fvals0(:,i).*fvals0(:,j);
  end
end
%
mydiff = abs(fvals - reshape(fvals0',[Norb npbox nboxes]));

%
nhess = ndim*(ndim+1)/2;
%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
nd = size(interpolating_vectors,2);
phi_kl = zeros(nd,npbox,nboxes); 
for k = 1:nd
  phi_kl(k,:,:) = reshape(interpolating_vectors(:,k),[npbox nboxes]);
end
pot = zeros(nd,npbox,nboxes);
nd0 = 10;
nchnk = floor(nd/nd0);
% for jchnk = 1:1
for jchnk = 1:nchnk
  jidxv = (jchnk-1)*nd0 + (1:nd0);
  phi_kl0 = phi_kl(jidxv,:,:);
  pot0=zeros(nd0,npbox,nboxes);
  grad0=zeros(nd0,ndim,npbox,nboxes);
  hess0=zeros(nd0,nhess,npbox,nboxes);
  tic
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
end
jchnk = nchnk+1;
jidxv = ((jchnk-1)*nd0+1):nd;
nd0 = numel(jidxv);
phi_kl0 = phi_kl(jidxv,:,:);
pot0=zeros(nd0,npbox,nboxes);
grad0=zeros(nd0,ndim,npbox,nboxes);
hess0=zeros(nd0,nhess,npbox,nboxes);
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

%
phi_ij_leaf = zeros(nd,npbox,nleafbox); % leaf box information
potleaf = zeros(nd,npbox,nleafbox);
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
    end
  end
end

%% compute smooth volume integral
% nd = nd0;
Vmunu = zeros(nd,nd);
% ell = 1;
for ell = 1:nd
  % (2*L)^2 * potleaf = (2*L)^2 * \int_{-1/2}^{1/2} (phi_k*phi_l)/|r' - r| dV
  potleaf_ell = potleaf(ell,:,:)/ratio^2; 
  for j = 1:nd
    % (2*L)^3 * \int_{-1/2}^{1/2} phi_{ij} * pot_{kl} dV
    Vmunu(j,ell) = sum(squeeze(phi_ij_leaf(j,:,:).*potleaf_ell).*wtsleaf(:,:),'all')/ratio^3;
  end
end


Vijkl = zeros(Norb,Norb,Norb,Norb);
tmpidx = 0;
for i = 1:Norb
  for j = 1:Norb
    for k = 1:Norb
      for l = 1:Norb
        for mu = 1:nd
          for nu = 1:nd
            Vijkl(i,j,k,l) = Vijkl(i,j,k,l) ...
                        + Vmunu(mu,nu)*collocation_matrix(i,mu)*collocation_matrix(j,mu) ...
                        *collocation_matrix(k,nu)*collocation_matrix(l,nu);
          end
        end
      end
    end
  end
end

save('Vmunu_h2o_ccpvdz.mat','Vmunu')
if exist("Vmunu_h2o_ccpvdz.h5", 'file') ~= 2
  h5create("Vmunu_h2o_ccpvdz.h5","/DS1",[nd nd])
  h5write("Vmunu_h2o_ccpvdz.h5","/DS1",Vmunu)
  h5disp("Vmunu_h2o_ccpvdz.h5")
end
keyboard