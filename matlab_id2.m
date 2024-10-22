% dummy test 2
% 10/21/24 Hai

addpath('./utils/')
addpath('./treefun/')

clear all
order = 10;
eps = 1e-03; 

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
% this have 5 digits accuracy...
% Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
%                 ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
%                 nboxes,nlevels,ltree,itree,iptr,centers,boxsize);

%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
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

%%%
ndsk = 285;
[SK,RD,T] = id(A,ndsk); % somewhere between 200 and 300 is a good number for 10 digits...
%
Ask = A(:,SK);
idcoefs = zeros(ndsk,size(A,2));
idcoefs(:,SK) = eye(ndsk);
idcoefs(:,RD) = T;
Aid = Ask*idcoefs;
diff = abs(Aid - A);
figure(1),clf,
imagesc(log10(diff)); colorbar, caxis([-15 -10])

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

% 12 digits accuracy
potid = tensorprod(Ask,potsk,2,1);
diff = abs(pot - potid);

keyboard
