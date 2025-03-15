%
%
% 02/27/25 Hai

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

% setup
% isdf_base_path = '/mnt/home/cyeh/ceph/papers/isdf_adaptive/H2O_dimer/ccpvdz/isdf_adap/';
bdmk_exec = '../../../utils/f/int2-bdmk-mlscf';
treefun_order = 8;
treefun_eps = 1e-06; 
isdf_eps = 1e-3;
nd = 2800;

% treefun_order = 6;
% treefun_eps = 1e-02; 

%%% resolve tree on cgto
rad = 15;
geom = sprintf([ ...
        'O  -1.551007  -0.114520   0.000000.\n',...
        'H  -1.934259   0.762503   0.000000.\n',...
        'H  -0.599677   0.040712   0.000000.\n',...
        'O   1.350625   0.111469   0.000000.\n',...
        'H   1.680398  -0.373741  -0.758561.\n',...
        'H   1.680398  -0.373741   0.758561.\n']),
molname = 'h2o_dimer';
basmod = 'cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
opts = struct('balance',true,...
              'tol',treefun_eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
% func = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
disp("=========Start treefun=======");
tic
f = treefun3(func,[-rad rad -rad rad -rad rad],treefun_order,opts);
time = toc;
disp("    treefun time is : " + time + " seconds");
disp("    treefun order is : " + treefun_order);
disp("    treefun num of boxes is : " + size(f.domain,2));
disp("=========End treefun=======");
% figure(1),clf,plot(f,func);

%%% treefun to bdmk
disp("=========Start upsample adaptive grid=======");
tic
Norb = mol.nao_nr; % 
ndim = 3;
ratio = 0.5/rad; % from boxlen to 1
ipoly = 0;
f2 = f; 
f2.n = floor(1.5*f.n);
[src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f2,ndim,ratio,ipoly);
if 1
  DS1 = src/ratio;
else
  DS1 = srcleaf/ratio;
end
time = toc;
disp("    upsample adaptive grid time is : " + time + " seconds");
disp("=========End upsample adaptive grid=======");
disp("    ");

%%% eval cgto
disp("=========Start Norb^2 basis eval=======");
tic
src0 = src/ratio;
npts = numel(src0(1,:));
fvals0 = squeeze(mol.eval_gto(eval_name, cat(4,squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:)))));
fvals0 = reshape(fvals0,[npts Norb]);
fvals_ij = zeros(npts,Norb*(Norb+1)/2);
tmpidx = 0;
for i=1:Norb
  for j=1:i
    tmpidx = tmpidx + 1;
    fvals_ij(:,tmpidx) = fvals0(:,i).*fvals0(:,j);
  end
end
time = toc;
disp("    Norb^2 basis eval time is : " + time + " seconds");
disp("=========End Norb^2 basis eval=======");
disp("    ");

%%% load isdf data
if 0
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
end

%%%%%%
disp("=========Start interpolative decomposition=======");
tic
A = fvals_ij;
tic
[SK,RD,T] = id(A,nd); %
Ask = A(:,SK);
idcoefs = zeros(nd,size(A,2));
idcoefs(:,SK) = eye(nd);
idcoefs(:,RD) = T;
% Aid = Ask*idcoefs;
toc
dims_rd = [2,npts,nd];
interpolating_vectors = zeros(dims_rd);
interpolating_vectors(1,:,:) = Ask;
eps_string = regexprep(sprintf('%.0e', isdf_eps), 'e-0*', 'e-'); 
isdf_filename = sprintf('myisdf_%s.h5', eps_string);
file_id = H5F.create(isdf_filename, "H5F_ACC_TRUNC", "H5P_DEFAULT", "H5P_DEFAULT"); 
% save interpolating_vectors
ipv_dataspace_id = H5S.create_simple(3, fliplr(size(interpolating_vectors)), []);
ipv_dcpl_id = H5P.create('H5P_DATASET_CREATE');
ipv_dataset_id = H5D.create(file_id, '/interpolating_vectors', 'H5T_NATIVE_DOUBLE', ipv_dataspace_id, ipv_dcpl_id);
H5D.write(ipv_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', interpolating_vectors);
H5D.close(ipv_dataset_id);
H5F.close(file_id);
time = toc;
disp("    id time is : " + time + " seconds");
disp("=========End interpolative decomposition=======");
disp("    ");

%%% write treefun data (bdmk format) to .h5
disp("=========Start bdmk=======");
tic
ikernel = 1;
beta = 6.0d0;
write_treefun_h5(molname,basmod,ndim,treefun_eps,ikernel,beta,...
                          ipoly,norder,npbox,nboxes,nlevels,ltree,itree,...
                          iptr,centers,boxsize,DS1,nleafbox,wtsleaf,ratio)
%%% call bdmk
treefun_filename = sprintf('treefun_%s_%s.h5', molname, erase(basmod, '.dat'));
output_filename = sprintf('bdmk_%s.h5', eps_string);
[status1, cmdout1] = system(sprintf('chmod +x %s', bdmk_exec));
[status2, cmdout2] = system(sprintf('%s %s %s %s', ...
                                    bdmk_exec, treefun_filename, isdf_filename, output_filename));
%
log_filename = 'bdmk_output.txt';  
fid = fopen(log_filename, 'w');
fprintf(fid, '%s\n', cmdout2);
fclose(fid);
fprintf('Saved bdmk output to: %s\n', log_filename);
time = toc;
disp("    bdmk time is : " + time + " seconds");
disp("=========End bdmk=======");
disp("    ");

%%% Vmunu
Vmunu = h5read(output_filename, '/Vmunu');

%%% Vijkl, after loading isdf_file and output_file
disp("=========Start Vijkl computation=======");
tic
Vijkl = zeros(Norb,Norb,Norb,Norb);
Vijkl = computeVijkl_mex(nd, Norb, idcoefs, Vmunu, Vijkl);
% Vijkl = computeVijkl(nd, Norb, idcoefs, Vmunu, Vijkl);
time = toc;
disp("    Vijkl time is : " + time + " seconds");
disp("=========End Vijkl computation=======");

% save
eri_mat_filename = ['ERI_' molname '_' erase(erase(basmod, '.dat'),'-') '_' eps_string '.mat'];
save(eri_mat_filename,'Vijkl')   
%
eri_h5_filename = ['ERI_' molname '_' erase(erase(basmod, '.dat'),'-') '_' eps_string '.h5'];
file_id = H5F.create(eri_h5_filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
vijkl_dims = [Norb, Norb, Norb, Norb];
vijkl_dataspace_id = H5S.create_simple(4, fliplr(vijkl_dims), []);
vijkl_dcpl_id = H5P.create('H5P_DATASET_CREATE');
vijkl_dataset_id = H5D.create(file_id, '/DS1', 'H5T_NATIVE_DOUBLE', vijkl_dataspace_id, vijkl_dcpl_id);
H5D.write(vijkl_dataset_id, 'H5T_NATIVE_DOUBLE', 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', Vijkl);
H5D.close(vijkl_dataset_id);
H5S.close(vijkl_dataspace_id);
H5F.close(file_id);
system(sprintf('rm -f %s', isdf_filename))

% keyboard
