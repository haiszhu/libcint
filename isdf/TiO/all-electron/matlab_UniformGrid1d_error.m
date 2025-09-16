%
%
% 09/16/25 Hai

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

% setup
% isdf_base_path = '/mnt/home/cyeh/ceph/papers/isdf_adaptive/H2O_dimer/ccpvdz/isdf_adap/';
bdmk_exec = '../../../utils/f/int2-bdmk-mlscf';
treefun_order = 4;
treefun_eps = 1e-03; 
% treefun_order = 7;
% treefun_eps = 1e-06; 
isdf_eps = 1e-3;
nd = 290;

%%% in the case of 'cc-pvdz-ccECP.dat', we need larger rad
%%% define mol
rad = 30;
geom = sprintf([ ...
    'O  0.000000  0.000000  -1.623000\n',...
    'Ti 0.000000  0.000000   0.000000\n']),
molname = 'TiO';
% basmod = 'cc-pvdz.dat';
basmod = 'cc-pvtz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
opts = struct('balance',true,...
              'tol',treefun_eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
func  = @(x,y,z)   mol.eval_gto( eval_name, cat(4,x,y,z));
funci = @(x,y,z,i) mol.eval_gtoi(eval_name, cat(4,x,y,z), i);
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
disp("=========Adaptive grid info=======");
disp("    User prescribed tolerance is : " + sprintf('%.3e', treefun_eps));
disp("    Adaptive grid size would be : " + numel(leaves(f))*(ceil(1.5*treefun_order))^3);
disp("    ");


disp("=========Verify uniform grid error=======");
disp("    This is actually a 1d slice... ");
Nall = 2.^(6:17);
for j=1:numel(Nall)
  %
  N = Nall(j);
  %
  [x] = ndgrid( (0:N-1)/N );      % unit cube [0,1)^3
  x = (x-1/2)*2*rad;
  y = 0*x;
  z = 0*x;
  %
  M = 2*N;
  [xm] = ndgrid( (0:M-1)/M );
  xm = (xm-1/2)*2*rad;
  ym = 0*xm;
  zm = 0*xm;
  w0m = (1-0)/M*ones(M,1)*2*rad;
  %
  Rel_Err = zeros(mol.nao_nr,1);
  for i = 1:mol.nao_nr
    %
    GTOival = funci(x,y,z,i);
    FN = fftn(GTOival);
    idx = (M/2-N/2+1):(M/2+N/2);
    FM = zeros(M,1,'like',FN);
    FM(idx) = fftshift(FN);
    FM = ifftshift(FM);
    %
    GTOivalM = real( ifft(FM,[M]) * (M/N) );
    %
    GTOivalM_ref = funci(xm,ym,zm,i);
    % % 
    % diff = abs(GTOivalM - GTOivalM_ref);
    % % 
    % Rel_Err(i) = max(diff(:))/max(max(abs(GTOivalM_ref(:))),1);
    %
    int2_Err = sqrt(sum((GTOivalM - GTOivalM_ref).^2.*w0m,'all'));
    int2     = sqrt(sum((GTOivalM_ref).^2.*w0m,'all'));
    % 
    Rel_Err(i) = int2_Err/max(int2,1);
    %
    % disp("    " + i + "-th basis eval error is : " + Rel_Err(i));
  end
  % disp("=========End of verifying uniform grid error=======");
  % disp("    ");
  disp("    Max relative L2 norm error is : " + max(Rel_Err));
  disp("    Uniform grid size would be " + sprintf('%d×%d×%d ', N, N, N));
  disp("    ");
end

keyboard

