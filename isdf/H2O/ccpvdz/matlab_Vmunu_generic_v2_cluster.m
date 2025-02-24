% h2o
% 
% module load hdf5
% module load matlab
%
% 02/24/25 Hai

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')
addpath('/mnt/home/cyeh/ceph/papers/isdf_adaptive/isdf/h2o/ccpvdz_src_1e-8')

clear all

isdf_base_path = '/mnt/home/cyeh/ceph/papers/isdf_adaptive/isdf/h2o/ccpvdz_src_1e-8';
bdmk_exec = '../../../utils/f/int2-bdmk-mlscf';
treefun_order = 10;
treefun_eps = 1e-8; 
isdf_eps = 1e-3;

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
func2 = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
f = treefun3(func2,[-rad rad -rad rad -rad rad],treefun_order,opts);
% plot(f,func2);

%%% treefun to bdmk
Norb = mol.nao_nr; % 
ndim = 3;
ratio = 0.5/rad; % from boxlen to 1
ipoly = 0;
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);
if 1
  DS1 = src/ratio;
else
  DS1 = srcleaf/ratio;
end

%%% write treefun data (bdmk format) to .h5
ikernel = 1;
beta = 6.0d0;
write_treefun_h5(molname,basmod,ndim,treefun_eps,ikernel,beta,...
                          ipoly,norder,npbox,nboxes,nlevels,ltree,itree,...
                          iptr,centers,boxsize,DS1,nleafbox,wtsleaf,ratio)

%%% call bdmk
eps_string = regexprep(sprintf('%.0e', isdf_eps), 'e-0*', 'e-'); 
treefun_filename = sprintf('treefun_%s_%s.h5', molname, erase(basmod, '.dat'));
% isdf_filename = sprintf('isdf_%s.h5', eps_string);
isdf_filename = sprintf('%s/isdf_%s.h5', isdf_base_path, eps_string);
output_filename = sprintf('bdmk_%s.h5', eps_string);
[status1, cmdout1] = system(sprintf('chmod +x %s', bdmk_exec));
[status2, cmdout2] = system(sprintf('%s %s %s %s', ...
                                    bdmk_exec, treefun_filename, isdf_filename, output_filename));
disp(cmdout2);
info = h5info(output_filename);
h5disp(output_filename)

%%% Vmunu
Vmunu = h5read(output_filename, '/Vmunu');

%%% Vijkl, after loading isdf_file and output_file
collocation_matrix = h5read(isdf_filename, '/collocation_matrix'); % this is not correct
collocation_matrix = squeeze(collocation_matrix(1,:,:));
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

% save
eri_mat_filename = ['ERI_h2o_ccpvdz_eps_' eps_string '.mat'];
save(eri_mat_filename,'Vijkl')    
%
eri_h5_filename = ['ERI_h2o_ccpvdz_eps_' eps_string '.h5'];
file_id = H5F.create(eri_h5_filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
vijkl_dims = [Norb, Norb, Norb, Norb];
vijkl_dataspace_id = H5S.create_simple(4, fliplr(vijkl_dims), []);
vijkl_dcpl_id = H5P.create('H5P_DATASET_CREATE');
vijkl_dataset_id = H5D.create(file_id, '/DS1', 'H5T_NATIVE_DOUBLE', vijkl_dataspace_id, vijkl_dcpl_id);
H5D.write(vijkl_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', Vijkl);
H5D.close(vijkl_dataset_id);
H5S.close(vijkl_dataspace_id);
H5F.close(file_id);
% %
% h5create(eri_h5_filename,"/DS1",[Norb Norb Norb Norb])
% h5write(eri_h5_filename,"/DS1",Vijkl)
% h5disp(eri_h5_filename)

% keyboard
