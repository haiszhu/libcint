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
order = 5;
eps = 1e-03; % 
% order = 8;
% eps = 1e-06;

%% geom
geom = sprintf([ ...
    'N  -1.578718  -0.046611   0.000000\n',...
    'H  -2.158621   0.136396  -0.809565\n',...
    'H  -2.158621   0.136396   0.809565\n',...
    'H  -0.849471   0.658193   0.000000\n',...
    'N   1.578718   0.046611   0.000000\n',...
    'H   2.158621  -0.136396  -0.809565\n',...
    'H   0.849471  -0.658193   0.000000\n',...
    'H   2.158621  -0.136396   0.809565\n']),

%% mol
basmod = 'cc-pvdz.dat';
% basmod = 'aug-cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);

%% eval
eval_name = 'GTOval_sph';
if exist('nh3_dimer_basis_check.mat','file')
  [x, y, z] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
  vals = mol.eval_gto(eval_name, cat(4,x,y,z));
  valsrs = reshape(vals,[],mol.nao_nr);
  load('nh3_dimer_basis_check.mat')
  diff = abs(vals - valsrs);
  max(diff(:))
end

%% adap tree
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
% func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
% f = treefun3(func,[-15 15 -15 15 -15 15],order,opts);
% plot(f,func);
func = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
f = treefun3(func,[-15 15 -15 15 -15 15],order,opts);
plot(f,func);

profile viewer

keyboard

