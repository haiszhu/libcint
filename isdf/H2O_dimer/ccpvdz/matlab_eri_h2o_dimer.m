%
%
% 01/20/25 Hai

% profile clear
% profile on

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all
order = 7;
eps = 1e-03; % 
% eps = 1e-04; %
% eps = 1e-05;

% create some verification mesh, supposedly to be the same as python_eri_h2o_dimer.py
% run python_eri_h2o_dimer.py to load x,y,z coordinates and basis from pyscf
if 1
  load('h2o_dimer_basis_check.mat')
  [x2, y2, z2] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
  xyz2 = [x2(:) y2(:) z2(:)];
  func = @(x,y,z) cgtofunc_h2o_dimer_ccpvdz(x,y,z);
  vals2 = func(x2, y2, z2);
  vals2rs = reshape(vals2,[],48);
  diff = abs(vals - reshape(vals2rs,[],48));
  max(diff(:))
end

% keyboard

%%% resolve tree on cgto^2
func2 = @(x,y,z) cgto2func_h2o_dimer_ccpvdz(x,y,z);
% func2 = @(x,y,z) cgtofunc_h2o_dimer_ccpvdz(x,y,z);
checkpts = [
          -1.551007  -0.114520   0.000000;... 
          -1.934259   0.762503   0.000000;... 
          -0.599677   0.040712   0.000000;... 
           1.350625   0.111469   0.000000;... 
           1.680398  -0.373741  -0.758561;... 
           1.680398  -0.373741   0.758561;... 
           ]';

opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',checkpts, ... 
              'ifcoeffs',false);
f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 
plot(f,func2);

%%% treefun to bdmk
Norb = 48;
ndim = 3;
ratio = 0.5/15; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

%%% eval cgto
src0 = src/ratio;
func = @(x,y,z) cgtofunc_h2o_dimer_ccpvdz(x,y,z);
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

%%% save data
eps_str = sprintf('%.0e', eps);
mat_filename = sprintf('ERI_h2o_dimer_ccpvdz_%s.mat', eps_str);
h5_filename = sprintf('ERI_h2o_dimer_ccpvdz_%s.h5', eps_str);
save(mat_filename,'Vijkl')    
if exist(h5_filename, 'file') ~= 2
  h5create(h5_filename,"/DS1",[Norb Norb Norb Norb])
end
h5write(h5_filename,"/DS1",Vijkl)
h5disp(h5_filename)

% profile viewer

keyboard
