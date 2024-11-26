%
%
%

addpath('./utils/')
addpath('./treefun/')

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
src = src/ratio;

h5create("src_h2o_ccpvdz.h5","/DS1",[3 npbox nboxes])
h5write("src_h2o_ccpvdz.h5","/DS1",src)
h5disp("src_h2o_ccpvdz.h5")

keyboard
