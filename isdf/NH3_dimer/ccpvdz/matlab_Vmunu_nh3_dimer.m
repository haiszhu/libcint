%
%
% 02/10/25 Hai

profile clear
profile on

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

%% specify eps and order
order = 4;
eps = 1e-02; % 
% order = 8;
% eps = 1e-06;

rad = 15;

%% geom
% geom = sprintf([ ...
%     'N  -1.578718  -0.046611   0.000000\n',...
%     'H  -2.158621   0.136396  -0.809565\n',...
%     'H  -2.158621   0.136396   0.809565\n',...
%     'H  -0.849471   0.658193   0.000000\n',...
%     'N   1.578718   0.046611   0.000000\n',...
%     'H   2.158621  -0.136396  -0.809565\n',...
%     'H   0.849471  -0.658193   0.000000\n',...
%     'H   2.158621  -0.136396   0.809565\n']),
geom = sprintf([ ...
    'O    0    0.       0.\n',...
    'H    0    -0.757   0.587\n',...
    'H    0    0.757    0.587\n']),

%% mol
basmod = 'cc-pvdz.dat';
% basmod = 'aug-cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);

%% eval
eval_name = 'GTOval_sph';
% if exist('nh3_dimer_basis_check.mat','file')
%   [x, y, z] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
%   vals = mol.eval_gto(eval_name, cat(4,x,y,z));
%   valsrs = reshape(vals,[],mol.nao_nr);
%   load('nh3_dimer_basis_check.mat')
%   diff = abs(vals - valsrs);
%   max(diff(:))
% end

%% adap tree
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
% func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
% f = treefun3(func,[-rad rad -rad rad -rad rad],order,opts);
% plot(f,func);
func = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
f = treefun3(func,[-rad rad -rad rad -rad rad],order,opts);
plot(f,func);

profile viewer

%% treefun to bdmk
Norb = mol.nao_nr;
ndim = 3;
ratio = 0.5/rad; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

%% eval cgto
src0 = src/ratio;
func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
fvals0 = squeeze(func(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
fvals0 = permute(fvals0,[3 1 2]);
fvals = fvals0;

%% save src to .h5 'DS1'
eps_str = sprintf('%.0e', eps);
h5_filename = sprintf('src_nh3_dimer_ccpvdz_%s.h5', eps_str);
%
src_save = false;
if src_save
  if exist(h5_filename, 'file') ~= 2
    h5create(h5_filename,"/DS1",[3 npbox nboxes])
  end
  h5write(h5_filename,"/DS1",src)
  h5disp(h5_filename)
end
%
src_load = false;
if src_load
  info = h5info('src_nh3_dimer_ccpvdz_1e-06.h5');
  src_data = h5read('src_nh3_dimer_ccpvdz_1e-06.h5', ['/' 'DS1']);
  diff = abs(src_data - src);
  max(diff(:))
end


%% compute V_ijkl, should do a "mol_nh3_dimer.intor('int2e')" like interface
nd = Norb*(Norb+1)/2;
ikernel = 1;
beta = 6.0d0;
% Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
%                 ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
%                 nboxes,nlevels,ltree,itree,iptr,centers,boxsize);


%% compute V_munu
nd = Norb^2;
phi_kl = zeros(nd,npbox,nboxes); 
idx = 0;
for k = 1:Norb
  for ell = 1:Norb
    idx = idx + 1;
    phi_kl(idx,:,:) = fvals(k,:,:).*fvals(ell,:,:);
  end
end
A = reshape(permute(phi_kl,[2 3 1]),[],Norb^2)';
ndsk = 200;
[SK,RD,T] = id(A,ndsk); % somewhere between 200 and 300 is a good number for 10 digits...
%
Ask = A(:,SK);
idcoefs = zeros(ndsk,size(A,2));
idcoefs(:,SK) = eye(ndsk);
idcoefs(:,RD) = T;
Aid = Ask*idcoefs;
diff = abs(Aid - A);
figure(1),clf,
imagesc(log10(diff)); colorbar, caxis([-15 log10(eps)])

%%% only SK volume potential
fvalssk = reshape(idcoefs,[ndsk npbox nboxes]);
ikernel = 1;
beta = 6.0d0;
nhess = ndim*(ndim+1)/2;
%
potsk=zeros(ndsk,npbox,nboxes);
gradsk=zeros(ndsk,ndim,npbox,nboxes);
hesssk=zeros(ndsk,nhess,npbox,nboxes);
tic
ifpgh=1;
ifpghtarg=0;
ntarg = 100;
targs=zeros(ndim,ntarg); pote=zeros(ndsk,ntarg);
grade=zeros(ndsk,ndim,ntarg); hesse=zeros(ndsk,nhess,ntarg);
timeinfo = zeros(20,1);
[potsk,gradsk,hesssk,pote,grade,hesse] = ...
bdmk_mex(ndsk,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvalssk, ...
       ifpgh,potsk,gradsk,hesssk,ntarg,targs, ...
       ifpghtarg,pote,grade,hesse,timeinfo);
potid = tensorprod(Ask,potsk,2,1);

%%% everything volume potential
nd = Norb^2;
phi_kl;
ikernel = 1;
beta = 6.0d0;
nhess = ndim*(ndim+1)/2;
%
pot=zeros(nd,npbox,nboxes);
grad=zeros(nd,ndim,npbox,nboxes);
hess=zeros(nd,nhess,npbox,nboxes);
tic
ifpgh=1;
ifpghtarg=0;
ntarg = 100;
targs=zeros(ndim,ntarg); pote=zeros(nd,ntarg);
grade=zeros(nd,ndim,ntarg); hesse=zeros(nd,nhess,ntarg);
timeinfo = zeros(20,1);
[pot,grad,hess,pote,grade,hesse] = ...
bdmk_mex(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,phi_kl, ...
       ifpgh,pot,grad,hess,ntarg,targs, ...
       ifpghtarg,pote,grade,hesse,timeinfo);

diff = abs(pot - potid);
max(diff(:))

figure(2),clf,
imagesc(log10(squeeze(diff(1,:,:)))); 
colorbar

keyboard

iname = 1;
idfnamei = './isdf_adap/isdf_1e-6.h5';
info = h5info(idfnamei);
Np = h5read(idfnamei, '/Np');
collocation_matrix = h5read(idfnamei, '/collocation_matrix');
interpolating_points = h5read(idfnamei, '/interpolating_points');
interpolating_vectors = h5read(idfnamei, '/interpolating_vectors');
kpts = h5read(idfnamei, '/kpts');
nkpts_ibz = h5read(idfnamei, '/nkpts_ibz');
nqpts_ibz = h5read(idfnamei, '/nqpts_ibz');
qpts = h5read(idfnamei, '/qpts');


% profile viewer

keyboard

