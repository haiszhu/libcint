%
%
% 02/11/25 Hai

profile clear
profile on

addpath('./utils/')
addpath('./treefun/')
addpath('./data')

clear all
order = 8;
eps = 1e-06; 

%%% resolve tree on cgto^2
rad = 15;
geom = sprintf([ ...
    'O    0    0.       0.\n',...
    'H    0    -0.757   0.587\n',...
    'H    0    0.757    0.587\n']),
basmod = 'cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), './basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
opts = struct('balance',true,...
              'tol',eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
func2 = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
f = treefun3(func2,[-rad rad -rad rad -rad rad],order,opts);
plot(f,func2);

%%% treefun to bdmk
Norb = 24; % 
ndim = 3;
ratio = 0.5/15; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

r = h5read('src_h2o_ccpvdz.h5', '/DS1');
npts = numel(r(1,:));
r2 = reshape(r,3,[]);
diff = abs(reshape(src/ratio,3,[]) - r2);

%%% eval cgto
src0 = src/ratio;
func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
fvals0 = squeeze(func(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
fvals0 = permute(fvals0,[3 1 2]);
fvals = fvals0;

%%% compute V_ijkl
nd = Norb*(Norb+1)/2;
ikernel = 1;
beta = 6.0d0;
% Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
%                 ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
%                 nboxes,nlevels,ltree,itree,iptr,centers,boxsize);

%%% load isdf data
%
idfname = 'isdf_1e-9.h5'; % 0.0026
idfname = 'isdf_1e-3.h5'; % 0.4476
info = h5info(idfname);
Np = h5read(idfname, '/Np');
collocation_matrix = h5read(idfname, '/collocation_matrix');
interpolating_points = h5read(idfname, '/interpolating_points');
interpolating_vectors = h5read(idfname, '/interpolating_vectors');
kpts = h5read(idfname, '/kpts');
nkpts_ibz = h5read(idfname, '/nkpts_ibz');
nqpts_ibz = h5read(idfname, '/nqpts_ibz');
qpts = h5read(idfname, '/qpts');
collocation_matrix = squeeze(collocation_matrix(1,:,:));
interpolating_vectors = squeeze(interpolating_vectors(1,:,:));

%%% compute
nd = size(interpolating_vectors,2);
Vmunu = Vmunucomp(Norb,ratio,interpolating_vectors,nleafbox,srcleaf,wtsleaf,...
                ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                nboxes,nlevels,ltree,itree,iptr,centers,boxsize);

%%%
Vijkl = zeros(Norb,Norb,Norb,Norb);
for i = 1:Norb
  for j = 1:Norb
    collocation_matrix_ij = collocation_matrix(i,:).*collocation_matrix(j,:);
    for k = 1:Norb
      for l = 1:Norb
        collocation_matrix_kl = collocation_matrix(k,:).*collocation_matrix(l,:);
        Vijkl(i,j,k,l) = sum(Vmunu.*( collocation_matrix_ij(:) ...
                                     .*collocation_matrix_kl(:)'),'all');
      end
    end
  end
end

save('Vmunu_h2o_ccpvdz.mat','Vmunu')
if exist("Vmunu_h2o_ccpvdz.h5", 'file') ~= 2
  h5create("Vmunu_h2o_ccpvdz.h5","/DS1",[nd nd])
  h5write("Vmunu_h2o_ccpvdz.h5","/DS1",Vmunu)
  h5disp("Vmunu_h2o_ccpvdz.h5")
end

save('ERI_h2o_ccpvdz.mat','Vijkl')    
if exist("ERI_h2o_ccpvdz.h5", 'file') ~= 2
  h5create("ERI_h2o_ccpvdz.h5","/DS1",[Norb Norb Norb Norb])
  h5write("ERI_h2o_ccpvdz.h5","/DS1",Vijkl)
  h5disp("ERI_h2o_ccpvdz.h5")
end

profile viewer

keyboard
