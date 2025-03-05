%
%
% 02/27/25 Hai

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

% setup
% isdf_base_path = '/mnt/home/cyeh/ceph/papers/isdf_adaptive/H2O_dimer/ccpvdz/isdf_adap/';
bdmk_exec = '../../../utils/f/int2-bdmk-mlscf';
treefun_order = 8;
treefun_eps = 1e-06; 
isdf_eps = 1e-3;
nd = 1600;

treefun_order = 6;
treefun_eps = 1e-02; 

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
f = treefun3(func,[-rad rad -rad rad -rad rad],treefun_order,opts);
disp("=========End treefun=======");
figure(1),clf,plot(f,func);

%%% treefun to bdmk
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

%%% eval cgto
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

%%%%%%
A = fvals_ij;
tic
[SK,RD,T] = id(A,nd); %
Ask = A(:,SK);
idcoefs = zeros(nd,size(A,2));
idcoefs(:,SK) = eye(nd);
idcoefs(:,RD) = T;
% Aid = Ask*idcoefs;

keyboard
