% generate adaptive tree with treefun
%
% 07/28/24 Hai

addpath('./utils/')
addpath('./treefun/')

clear all
order = 7;
% eps = 1e-03; 
% eps = 1e-04;
eps = 1e-05;
% eps = 1e-06;
% eps = 1e-07;

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
Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
                ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                nboxes,nlevels,ltree,itree,iptr,centers,boxsize);

save('ERI_h2o_ccpvdz.mat','Vijkl')    
if exist("ERI_h2o_ccpvdz.h5", 'file') ~= 2
  h5create("ERI_h2o_ccpvdz.h5","/DS1",[Norb Norb Norb Norb])
  h5write("ERI_h2o_ccpvdz.h5","/DS1",Vijkl)
  h5disp("ERI_h2o_ccpvdz.h5")
end

keyboard