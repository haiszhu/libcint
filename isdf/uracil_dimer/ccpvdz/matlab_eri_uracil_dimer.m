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
  load('uracil_dimer_basis_check.mat')
  [x2, y2, z2] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
  xyz2 = [x2(:) y2(:) z2(:)];
  func = @(x,y,z) cgtofunc_uracil_dimer_ccpvdz(x,y,z);
  vals2 = func(x2, y2, z2);
  vals2rs = reshape(vals2,[],264);
  diff = abs(vals - reshape(vals2rs,[],264));
  max(diff(:))
end

% keyboard

%%% resolve tree on cgto^2
func2 = @(x,y,z) cgto2func_uracil_dimer_ccpvdz(x,y,z);
% func2 = @(x,y,z) cgtofunc_uracil_dimer_ccpvdz(x,y,z);
checkpts = [
        -1.4663316    1.0121693    0.0000000;... 
        -0.6281464    1.9142678    0.0000000;... 
         0.7205093    1.6882688    0.0000000;... 
         1.6367290    2.7052764    0.0000000;... 
         1.2769036    4.0061763    0.0000000;... 
        -0.1286005    4.3621549    0.0000000;... 
        -0.9777230    3.2396433    0.0000000;... 
        -0.5972229    5.4864066    0.0000000;... 
         2.0103504    4.7938642    0.0000000;... 
         1.0232515    0.7061820    0.0000000;... 
        -1.9700268    3.4323850    0.0000000;... 
         2.6690620    2.3883417    0.0000000;... 
         1.4663316   -1.0121693    0.0000000;... 
         0.6281464   -1.9142678    0.0000000;... 
        -0.7205093   -1.6882688    0.0000000;... 
        -1.6367290   -2.7052764    0.0000000;... 
        -1.2769036   -4.0061763    0.0000000;... 
         0.1286005   -4.3621549    0.0000000;... 
         0.9777230   -3.2396433    0.0000000;... 
         0.5972229   -5.4864066    0.0000000;... 
        -2.0103504   -4.7938642    0.0000000;... 
        -1.0232515   -0.7061820    0.0000000;... 
         1.9700268   -3.4323850    0.0000000;... 
        -2.6690620   -2.3883417    0.0000000];

opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',checkpts, ... 
              'ifcoeffs',false);
f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 
plot(f,func2);

%%% treefun to bdmk
Norb = 264;
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
