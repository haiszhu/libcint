%
%
% 01/20/25 Hai

profile clear
profile on

addpath('../')
addpath('../utils/')
addpath('../treefun/')

clear all
order = 7;
eps = 1e-03; % 
eps = 1e-05; % 

% create some verification mesh, supposedly to be the same as python_eri_nh3_dimer.py
% run python_eri_nh3_dimer.py to load x,y,z coordinates and basis from pyscf
if 1
  load('nh3_dimer_basis_check.mat')
  [x2, y2, z2] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
  xyz2 = [x2(:) y2(:) z2(:)];
  func = @(x,y,z) cgtofunc_nh3_dimer_ccpvdz(x,y,z);
  vals2 = func(x2, y2, z2);
  vals2rs = reshape(vals2,[],58);
  diff = abs(vals - reshape(vals2rs,[],58));
  max(diff(:))
end

% keyboard

%%% resolve tree on cgto^2
% func2 = @(x,y,z) cgto2func_h2o_ccpvtz(x,y,z);
func2 = @(x,y,z) cgtofunc_nh3_dimer_ccpvdz(x,y,z);
checkpts = [
          -1.578718  -0.046611   0.000000;... 
          -2.158621   0.136396  -0.809565;...
          -2.158621   0.136396   0.809565;...
          -0.849471   0.658193   0.000000;...
           1.578718   0.046611   0.000000;...
           2.158621  -0.136396  -0.809565;...
           0.849471  -0.658193   0.000000;...
           2.158621  -0.136396   0.809565
           ]';
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',checkpts, ... 
              'ifcoeffs',false);
f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 
plot(f,func2);

profile viewer

keyboard
