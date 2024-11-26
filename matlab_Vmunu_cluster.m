% 
% LD_PRELOAD=/mnt/sw/nix/store/kkhxdzyg5kaxrxrl6j7pwapq7v7xacc7-openblas-0.3.26/lib/libopenblas.so matlab
% 
% need to force matlab to use openblas...
%
% to submit a job
% module load openblas
% OPENBLAS_DIR=$(echo $PATH | tr ':' '\n' | grep "openblas")
% OPENBLAS_LIB="$OPENBLAS_DIR/../lib/libopenblas.so"
% LD_PRELOAD=$OPENBLAS_LIB matlab
%
% 11/26/24 Hai

addpath('./utils/')
addpath('./treefun/')
addpath('./data')
addpath('/mnt/home/cyeh/ceph/papers/isdf_adaptive/isdf/h2o/ccpvdz_src_1e-8')

clear all
order = 10;
eps = 1e-08; 

%%% resolve tree on cgto^2
func2 = @(x,y,z) cgto2func(x,y,z);
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',[0 0 0; 0 -0.757 0.757;0 0.587 0.587], ...
              'ifcoeffs',false);
f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 
% f = treefun3(func,[-0.5 0.5 -0.5 0.5 -0.5 0.5],order,opts); 
% plot(f,func2);

%%% treefun to bdmk
Norb = 24; % 
ndim = 3;
ratio = 0.5/15; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

r = h5read('src_h2o_ccpvdz_1e-08.h5', '/DS1');
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

idfname{1} = 'isdf_1e-3.h5'; % 0.0026
idfname{2} = 'isdf_1e-4.h5';
idfname{3} = 'isdf_1e-5.h5';
idfname{4} = 'isdf_1e-6.h5';
idfname{5} = 'isdf_1e-7.h5';
idfname{6} = 'isdf_1e-8.h5';
idfname{7} = 'isdf_1e-9.h5';
% idfname = 'isdf_1e-4.h5'; % 0.4476
for iname = 1:7
  idfnamei = idfname{iname};
  info = h5info(idfnamei);
  Np = h5read(idfnamei, '/Np');
  collocation_matrix = h5read(idfnamei, '/collocation_matrix');
  interpolating_points = h5read(idfnamei, '/interpolating_points');
  interpolating_vectors = h5read(idfnamei, '/interpolating_vectors');
  kpts = h5read(idfnamei, '/kpts');
  nkpts_ibz = h5read(idfnamei, '/nkpts_ibz');
  nqpts_ibz = h5read(idfnamei, '/nqpts_ibz');
  qpts = h5read(idfnamei, '/qpts');
    
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
  
  % extra eps info in filename
  [~, idbasename, ~] = fileparts(idfnamei);
  eps_string = extractAfter(idbasename, '_'); 
  
  % save Vmunu to matlab .mat
  mat_filename = ['Vmunu_h2o_ccpvdz_eps_' eps_string '.mat'];
  save(mat_filename, 'Vmunu');
  
  % save Vmunu to .h5
  h5_filename = ['Vmunu_h2o_ccpvdz_eps_' eps_string '.h5'];
  if exist(h5_filename, 'file') ~= 2
    h5create(h5_filename,"/DS1",[nd nd])
    h5write(h5_filename,"/DS1",Vmunu)
    h5disp(h5_filename)
  end
  
  % save Vijkl to matlab .mat
  eri_mat_filename = ['ERI_h2o_ccpvdz_eps_' eps_string '.mat'];
  save(eri_mat_filename,'Vijkl') 
  
  % save Vijkl to .h5
  eri_h5_filename = ['ERI_h2o_ccpvdz_eps_' eps_string '.h5'];
  if exist(eri_h5_filename, 'file') ~= 2
    h5create(eri_h5_filename,"/DS1",[Norb Norb Norb Norb])
    h5write(eri_h5_filename,"/DS1",Vijkl)
    h5disp(eri_h5_filename)
  end

end

% keyboard