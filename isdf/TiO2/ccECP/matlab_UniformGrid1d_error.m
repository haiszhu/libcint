%
%
% 07/22/25 Hai

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

% setup
% isdf_base_path = '/mnt/home/cyeh/ceph/papers/isdf_adaptive/H2O_dimer/ccpvdz/isdf_adap/';
bdmk_exec = '../../../utils/f/int2-bdmk-mlscf';
treefun_order = 8;
treefun_eps = 1e-04; 
isdf_eps = 1e-3;
nd = 290;

%%% define mol
rad = 15;
geom = sprintf([ ...
    'Ti    0    0.       0.\n',...
    'O     1.6  0.       0.\n',...
    'O    -1.6  0.       0.\n']),
molname = 'TiO2';
% basmod = 'cc-pvdz.dat';
basmod = 'cc-pvdz-ccECP.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
opts = struct('balance',true,...
              'tol',treefun_eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
func  = @(x,y,z)   mol.eval_gto( eval_name, cat(4,x,y,z));
% func = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
disp("=========Start treefun=======");
tic
f = treefun3(func,[-rad rad -rad rad -rad rad],treefun_order,opts);
time = toc;
disp("    treefun time is : " + time + " seconds");
disp("    treefun order is : " + treefun_order);
disp("    treefun num of boxes is : " + size(f.domain,2));
disp("=========End treefun=======");
disp("    ");
% figure(1),clf,plot(f,func);

keyboard

