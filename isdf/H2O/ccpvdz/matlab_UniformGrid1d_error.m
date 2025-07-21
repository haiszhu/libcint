% ICOSAHOM...
%
% 07/17/25

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

%%% define mol
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
func  = @(x,y,z)   mol.eval_gto( eval_name, cat(4,x,y,z));
funci = @(x,y,z,i) mol.eval_gtoi(eval_name, cat(4,x,y,z), i);
%
N = 64;
N = 128;
N = 256;
N = 512;
N = 1024;
N = 2048;
N = 4096;
% N = 8192;
% N = 16384;
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
%
Rel_Err = zeros(mol.nao_nr,1);
disp("=========Verify uniform grid error=======");
disp("    With uniform grid size : " + N);
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
  % 
  diff = abs(GTOivalM - GTOivalM_ref);
  % 
  Rel_Err(i) = max(diff(:))/max(max(abs(GTOivalM_ref(:))),1);
  %
  disp("    " + i + "-th basis eval error is : " + Rel_Err(i));
end
disp("=========End of verifying uniform grid error=======");
disp("    ");
disp("    Max relative error is : " + max(Rel_Err));


keyboard
