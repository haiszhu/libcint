% h2o
%
% 02/16/25 Hai


addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all
order = 10;
eps = 1e-08; 

%%% resolve tree on cgto^2
rad = 15;
geom = sprintf([ ...
    'O    0    0.       0.\n',...
    'H    0    -0.757   0.587\n',...
    'H    0    0.757    0.587\n']),
basmod = 'cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
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
Norb = mol.nao_nr; % 
ndim = 3;
ratio = 0.5/rad; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

%
r = h5read('src_h2o_ccpvdz_1e-08.h5', '/DS1');
diff = abs(src/ratio - r);

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
idfname = 'isdf_1e-3.h5'; % 0.8859462695394982
idfname = 'isdf_1e-4.h5'; % 0.4503030075599318
idfname = 'isdf_1e-5.h5'; % 0.31833810509011884
idfname = 'isdf_1e-6.h5'; % 
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

% extract
[~, idbasename, ~] = fileparts(idfname);
eps_string = extractAfter(idbasename, '_'); 

% save
mat_filename = ['Vmunu_h2o_ccpvdz_eps_' eps_string '.mat'];
save(mat_filename, 'Vmunu');
h5_filename = ['Vmunu_h2o_ccpvdz_eps_' eps_string '.h5'];
h5create(h5_filename,"/DS1",[nd nd])
h5write(h5_filename,"/DS1",Vmunu)
h5disp(h5_filename)

% save
eri_mat_filename = ['ERI_h2o_ccpvdz_eps_' eps_string '.mat'];
save(eri_mat_filename,'Vijkl')    
eri_h5_filename = ['ERI_h2o_ccpvdz_eps_' eps_string '.h5'];
h5create(eri_h5_filename,"/DS1",[Norb Norb Norb Norb])
h5write(eri_h5_filename,"/DS1",Vijkl)
h5disp(eri_h5_filename)


keyboard

