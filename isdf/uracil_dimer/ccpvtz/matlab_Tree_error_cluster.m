% 
%
%

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

epsvals = [1e-2 1e-3 1e-4 1e-5 1e-6];
ordervals = [3 4 5 6 7];

%
rad = 15;
geom = sprintf([ ...
    'O    -1.4663316    1.0121693    0.0000000.\n',...
    'C    -0.6281464    1.9142678    0.0000000.\n',...
    'N     0.7205093    1.6882688    0.0000000.\n',...
    'C     1.6367290    2.7052764    0.0000000.\n',...
    'C     1.2769036    4.0061763    0.0000000.\n',...
    'C    -0.1286005    4.3621549    0.0000000.\n',...
    'N    -0.9777230    3.2396433    0.0000000.\n',...
    'O    -0.5972229    5.4864066    0.0000000.\n',...
    'H     2.0103504    4.7938642    0.0000000.\n',...
    'H     1.0232515    0.7061820    0.0000000.\n',...
    'H    -1.9700268    3.4323850    0.0000000.\n',...
    'H     2.6690620    2.3883417    0.0000000.\n',...
    'O     1.4663316   -1.0121693    0.0000000.\n',...
    'C     0.6281464   -1.9142678    0.0000000.\n',...
    'N    -0.7205093   -1.6882688    0.0000000.\n',...
    'C    -1.6367290   -2.7052764    0.0000000.\n',...
    'C    -1.2769036   -4.0061763    0.0000000.\n',...
    'C     0.1286005   -4.3621549    0.0000000.\n',...
    'N     0.9777230   -3.2396433    0.0000000.\n',...
    'O     0.5972229   -5.4864066    0.0000000.\n',...
    'H    -2.0103504   -4.7938642    0.0000000.\n',...
    'H    -1.0232515   -0.7061820    0.0000000.\n',...
    'H     1.9700268   -3.4323850    0.0000000.\n',...
    'H    -2.6690620   -2.3883417    0.0000000\n']),
molname = 'uracil_dimer';
basmod = 'cc-pvtz.dat';
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
  relerr = treedata_error(f,func);
  disp("========= Relative L infinity error on order p grid =======");
  disp("    Error on Norb basis : " + relerr );
  disp("    ");
  relerrs(i) = relerr;

  %=================== err on order p grid, Norb^2 basis ====================
  % f2 = treedata_resample(f,func2,treefun_order);
  % relerr2 = treedata_error(f2,func2);
  relerr2 = treedata_resample_error(f,func2,treefun_order);
  disp("========= Relative L infinity error on order p grid =======");
  disp("    Error on Norb^2 basis : " + relerr2 );
  disp("    ");
  relerr2s(i) = relerr2;
  
  %=================== err on order 2*p grid, Norb^2 basis ==================
  upfactor = 1.5;
  % f2up = treedata_resample(f,func2,ceil(upfactor*treefun_order));
  % relerr2up = treedata_error(f2up,func2);
  relerr2up = treedata_resample_error(f,func2,ceil(upfactor*treefun_order));
  disp("========= Relative L infinity error on order 2*p grid =======");
  disp("    Error on Norb^2 basis : " + relerr2up );
  disp("    ");
  relerr2ups(i) = relerr2up;

  % keyboard
end

[~, basname, ~] = fileparts(basmod);
mat_filename = [molname '_' basname '.mat'];
save(mat_filename, 'relerrs','relerr2s','relerr2ups','epsvals','ordervals');

% keyboard
