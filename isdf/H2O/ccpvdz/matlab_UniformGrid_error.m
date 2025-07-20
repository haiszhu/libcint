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
% N = 256;
% N = 384;
[x,y,z] = ndgrid( (0:N-1)/N );      % unit cube [0,1)^3
x = (x-1/2)*2*rad;
y = (y-1/2)*2*rad;
z = (z-1/2)*2*rad;
%
M = 2*N;
[xm,ym,zm] = ndgrid( (0:M-1)/M );
xm = (xm-1/2)*2*rad;
ym = (ym-1/2)*2*rad;
zm = (zm-1/2)*2*rad;
%
Rel_Err = zeros(mol.nao_nr,1);
disp("=========Verify uniform grid error=======");
disp("    With uniform grid size : " + N);
for i = 1:mol.nao_nr
  %
  GTOival = funci(x,y,z,i);
  FN = fftn(GTOival);
  idx = (M/2-N/2+1):(M/2+N/2);
  FM = zeros(M,M,M,'like',FN);
  FM(idx,idx,idx) = fftshift(FN);
  FM = ifftshift(FM);
  %
  GTOivalM = real( ifftn(FM,[M M M]) * (M/N)^3 );
  %
  GTOivalM_ref = funci(xm,ym,zm,i);
  % 
  diff = abs(GTOivalM - GTOivalM_ref);
  % 
  Rel_Err(i) = max(diff(:))/max(abs(GTOivalM_ref(:)));
  %
  disp("    " + i + "-th basis eval error is : " + Rel_Err(i));
end
disp("=========End of verifying uniform grid error=======");
disp("    ");
disp("    Max relative error is : " + max(Rel_Err));


keyboard
