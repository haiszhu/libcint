% 
%
%

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

treefun_order = 6;
treefun_eps = 1e-04; 
%%% resolve tree on cgto^2
rad = 15;
geom = sprintf([ ...
    'O    0    0.       0.\n',...
    'H    0    -0.757   0.587\n',...
    'H    0    0.757    0.587\n']),
molname = 'h2o';
basmod = 'cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
opts = struct('balance',true,...
              'tol',treefun_eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
func2 = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
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

%=================== err on order p grid, Norb basis ======================
relerr = treedata_error(f,func);
disp("========= Relative L infinity error on order p grid =======");
disp("    Error on Norb basis : " + relerr );
disp("    ");

%=================== err on order p grid, Norb^2 basis ====================
f2 = treedata_resample(f,func2,treefun_order);
relerr2 = treedata_error(f2,func2);
disp("========= Relative L infinity error on order p grid =======");
disp("    Error on Norb^2 basis : " + relerr2 );
disp("    ");

%=================== err on order 2*p grid, Norb^2 basis ==================
upfactor = 2;
f2up = treedata_resample(f,func2,upfactor*treefun_order);
relerr2up = treedata_error(f2up,func2);
disp("========= Relative L infinity error on order 2*p grid =======");
disp("    Error on Norb^2 basis : " + relerr2up );
disp("    ");

keyboard