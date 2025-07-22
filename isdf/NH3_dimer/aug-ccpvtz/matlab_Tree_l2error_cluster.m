% 
%
%

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

epsvals = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
ordervals = [3 4 5 6 7 8 9];

%
rad = 15;
geom = sprintf([ ...
    'N  -1.578718  -0.046611   0.000000.\n',...
    'H  -2.158621   0.136396  -0.809565.\n',...
    'H  -2.158621   0.136396   0.809565.\n',...
    'H  -0.849471   0.658193   0.000000.\n',...
    'N   1.578718   0.046611   0.000000.\n',...
    'H   2.158621  -0.136396  -0.809565.\n',...
    'H   0.849471  -0.658193   0.000000.\n',...
    'H   2.158621  -0.136396   0.809565\n']),
molname = 'nh3_dimer';
basmod = 'aug-cc-pvtz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
func2 = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));

%
relerrs = zeros(size(epsvals));
relerr2s = zeros(size(epsvals));
relerr2ups = zeros(size(epsvals));
for i = 1:numel(epsvals)
  %
  treefun_order = ordervals(i);
  treefun_eps = epsvals(i);
  %
  opts = struct('balance',true,...
                'tol',treefun_eps, ...
                'checkpts',mol.checkpts, ... 
                'ifcoeffs',false);
  %
  disp("=========Start treefun=======");
  tic
  f = treefun3(func,[-rad rad -rad rad -rad rad],treefun_order,opts);
  time = toc;
  disp("    treefun time is : " + time + " seconds");
  disp("    treefun order is : " + treefun_order);
  disp("    treefun num of boxes is : " + size(f.domain,2));
  disp("=========End treefun=======");
  disp("    ");

  %=================== err on order p grid, Norb basis ======================
  relerr = treedata_l2error(f,func);
  disp("========= Relative L^2 error on order p grid =======");
  disp("    Error on Norb basis : " + relerr );
  disp("    ");
  relerrs(i) = relerr;

  %=================== err on order p grid, Norb^2 basis ====================
  % f2 = treedata_resample(f,func2,treefun_order);
  % relerr2 = treedata_error(f2,func2);
  relerr2 = treedata_resample_l2error(f,func2,treefun_order);
  disp("========= Relative L^2 error on order p grid =======");
  disp("    Error on Norb^2 basis : " + relerr2 );
  disp("    ");
  relerr2s(i) = relerr2;
  
  %=================== err on order 2*p grid, Norb^2 basis ==================
  upfactor = 1.5;
  % f2up = treedata_resample(f,func2,ceil(upfactor*treefun_order));
  % relerr2up = treedata_error(f2up,func2);
  relerr2up = treedata_resample_l2error(f,func2,ceil(upfactor*treefun_order));
  disp("========= Relative L^2 error on order 2*p grid =======");
  disp("    Error on Norb^2 basis : " + relerr2up );
  disp("    ");
  relerr2ups(i) = relerr2up;

  % keyboard
end

[~, basname, ~] = fileparts(basmod);
mat_filename = [molname '_' basname '.mat'];
save(mat_filename, 'relerrs','relerr2s','relerr2ups','epsvals','ordervals');

% keyboard
