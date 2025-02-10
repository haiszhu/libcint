%
%
% 02/10/25 Hai

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

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
mol_nh3_dimer = gto(geom,basis);

%% eval
eval_name = 'GTOval_sph';
[x, y, z] = ndgrid(linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
vals = mol_nh3_dimer.eval_gto(eval_name, cat(4,x,y,z));
valsrs = reshape(vals,[],mol_nh3_dimer.nao_nr);
load('nh3_dimer_basis_check.mat')
diff = abs(vals - valsrs);
max(diff(:))

keyboard

% specify eps and order
% order = 5;
% eps = 1e-03; % 
order = 8;
eps = 1e-06;


keyboard

