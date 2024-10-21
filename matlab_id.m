% dummy id test
%
% 10/20/24 Hai 

addpath('./utils/')
addpath('./treefun/')

clear all

% learn id
n = 1000;
A = hilb(n);
rank(A)
[SK,RD,T] = id(A,20);
diff = abs(A(:,RD) - A(:,SK)*T);
figure(1),clf,
imagesc(log10(diff)); colorbar

% learn id 2
Norb = 40;
func3d = @(x,y,z,v) sin(x).*cos(y).*sin(z).*cos(v) ...
                    +cos(2*x).*sin(3*y).*cos(z).*sin(v);
[xx,yy,zz,vv] = ndgrid(linspace(-1,1,Norb));
xx = reshape(xx,[],Norb); yy = reshape(yy,[],Norb); zz = reshape(zz,[],Norb);
vv = reshape(vv,[],Norb);
val3d = func3d(xx,yy,zz,vv);
[SK,RD,T] = id(val3d,2);
diff = abs(val3d(:,RD) - val3d(:,SK)*T);
figure(1),clf,
imagesc(log10(diff)); colorbar % visualize id results...
val3dsk = val3d(:,SK); % useful stuff

% define sk func
vvsk = vv(1,SK);
func3dsk = @(x,y,z) tensorprod(sin(x).*cos(y).*sin(z),cos(vvsk)) ...
                    +tensorprod(cos(2*x).*sin(3*y).*cos(z),sin(vvsk));
eps = 1e-06;
order = 5;
opts = struct('balance',true,...
              'tol',eps, ...
              'ifcoeffs',false);
f = treefun3(func3dsk,[-1/2 1/2 -1/2 1/2 -1/2 1/2],order,opts); 
% f = treefun3(func,[-0.5 0.5 -0.5 0.5 -0.5 0.5],order,opts); 
plot(f,func3dsk);

%%% treefun to bdmk
ndim = 3;
ratio = 0.5/0.5; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);
src0 = src/ratio;

%%% skeleton computation
% eval sk func
fvalssk = squeeze(func3dsk(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
fvalssk = permute(fvalssk,[3 1 2]);

% compute V_sk of sk func
nd = numel(SK);
ikernel = 1;
beta = 6.0d0;
nhess = ndim*(ndim+1)/2;
%
potsk=zeros(nd,npbox,nboxes);
gradsk=zeros(nd,ndim,npbox,nboxes);
hesssk=zeros(nd,nhess,npbox,nboxes);
tic
ifpgh=1;
ifpghtarg=0;
ntarg = 100;
targs=zeros(ndim,ntarg); pote=zeros(nd,ntarg);
grade=zeros(nd,ndim,ntarg); hesse=zeros(nd,nhess,ntarg);
timeinfo = zeros(20,1);
[potsk,gradsk,hesssk,pote,grade,hesse] = ...
bdmk_mex(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvalssk, ...
       ifpgh,potsk,gradsk,hesssk,ntarg,targs, ...
       ifpghtarg,pote,grade,hesse,timeinfo);

%%% redundant computation
% use redundant stuff
ndrd = 15;
vvrd = vv(1,RD(1:ndrd));
func3drd = @(x,y,z) tensorprod(sin(x).*cos(y).*sin(z),cos(vvrd)) ...
                    +tensorprod(cos(2*x).*sin(3*y).*cos(z),sin(vvrd));

% eval rd func
fvalsrd = squeeze(func3drd(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
fvalsrd = permute(fvalsrd,[3 1 2]);

% compute V_rd of rd func
ndrd;
ikernel = 1;
beta = 6.0d0;
nhess = ndim*(ndim+1)/2;
%
potrd=zeros(ndrd,npbox,nboxes);
gradrd=zeros(ndrd,ndim,npbox,nboxes);
hessrd=zeros(ndrd,nhess,npbox,nboxes);
tic
ifpgh=1;
ifpghtarg=0;
ntarg = 100;
targs=zeros(ndim,ntarg); pote=zeros(ndrd,ntarg);
grade=zeros(ndrd,ndim,ntarg); hesse=zeros(ndrd,nhess,ntarg);
timeinfo = zeros(20,1);
[potrd,gradrd,hessrd,pote,grade,hesse] = ...
bdmk_mex(ndrd,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
       nboxes,nlevels,ltree,itree,iptr,centers,boxsize,fvalsrd, ...
       ifpgh,potrd,gradrd,hessrd,ntarg,targs, ...
       ifpghtarg,pote,grade,hesse,timeinfo);

%%% verify function value on tree
fvalserr = abs(permute(fvalsrd,[2 3 1]) - tensorprod(permute(fvalssk,[2 3 1]),T(:,1:ndrd),3,1));
fvalserr = reshape(fvalserr,[],ndrd);
figure(2),clf,
imagesc(log10(fvalserr)); colorbar
title('dummy func value err on grid: direct vs id')

%%% verify volume integral on tree
vvalserr = abs(permute(potrd,[2 3 1]) - tensorprod(permute(potsk,[2 3 1]),T(:,1:ndrd),3,1));
vvalserr = reshape(vvalserr,[],ndrd);
figure(3),clf,
imagesc(log10(vvalserr)); colorbar
title('volume integral value err on grid: direct vs id')

keyboard
