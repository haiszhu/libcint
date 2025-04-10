function f = cgto2func_h2o_dimer_aug_ccpvdz(x,y,z)
% implement cgto, which calls libcint, libcgto, np_helper, etc...
% phi_k * phi_l, totoal Norb^2, we store about half of it
% 
%
% h2o_dimer, aug-ccpvdz
% geom='''
%         O  -1.551007  -0.114520   0.000000
%         H  -1.934259   0.762503   0.000000
%         H  -0.599677   0.040712   0.000000
%         O   1.350625   0.111469   0.000000
%         H   1.680398  -0.373741  -0.758561
%         H   1.680398  -0.373741   0.758561
%         '''
% 
% mol_h2o_dimer = gto.M(
%         verbose=7, 
%         atom=geom, 
%         basis='aug-ccpvdz')
%
% use parameters from pyscf/pyscf/gto/eval_gto.py
% 
%
% 01/23/25 Hai

% respect the shape of input variables
Norb = 82;
nd = Norb*(Norb+1)/2;
[n1,n2,n3] = size(x);
f = zeros([n1 n2 n3 nd]); % need to be changed later

persistent ngrids shls_slice ao_loc non0tab atm bas env natm nbas nao triuidx triflag;
if isempty(ngrids)
  % ngrids = 1; 
  shls_slice = [0, 36];
  ao_loc = [ 0,  2,  3,  4,  7, 10, 13, 18, 23, 24, 25, 26, 29, 32, 33, 34, 35, 38, 41, 43, 44, 45, 48, 51, 54, 59, 64, 65, 66, 67, 70, 73, 74, 75, 76, 79, 82];
  non0tab = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
  atm = [
    8, 20,  1, 23,  0,  0,...
    1, 24,  1, 27,  0,  0,...
    1, 28,  1, 31,  0,  0,...
    8, 32,  1, 35,  0,  0,...
    1, 36,  1, 39,  0,  0,...
    1, 40,  1, 43,  0,  0];
  bas = [
    0,  0,  8,  2,  0, 58, 66,  0,...
    0,  0,  1,  1,  0, 82, 83,  0,...
    0,  0,  1,  1,  0, 84, 85,  0,...
    0,  1,  3,  1,  0, 86, 89,  0,...
    0,  1,  1,  1,  0, 92, 93,  0,...
    0,  1,  1,  1,  0, 94, 95,  0,...
    0,  2,  1,  1,  0, 96, 97,  0,...
    0,  2,  1,  1,  0, 98, 99,  0,...
    1,  0,  3,  1,  0, 44, 47,  0,...
    1,  0,  1,  1,  0, 50, 51,  0,...
    1,  0,  1,  1,  0, 52, 53,  0,...
    1,  1,  1,  1,  0, 54, 55,  0,...
    1,  1,  1,  1,  0, 56, 57,  0,...
    2,  0,  3,  1,  0, 44, 47,  0,...
    2,  0,  1,  1,  0, 50, 51,  0,...
    2,  0,  1,  1,  0, 52, 53,  0,...
    2,  1,  1,  1,  0, 54, 55,  0,...
    2,  1,  1,  1,  0, 56, 57,  0,...
    3,  0,  8,  2,  0, 58, 66,  0,...
    3,  0,  1,  1,  0, 82, 83,  0,...
    3,  0,  1,  1,  0, 84, 85,  0,...
    3,  1,  3,  1,  0, 86, 89,  0,...
    3,  1,  1,  1,  0, 92, 93,  0,...
    3,  1,  1,  1,  0, 94, 95,  0,...
    3,  2,  1,  1,  0, 96, 97,  0,...
    3,  2,  1,  1,  0, 98, 99,  0,...
    4,  0,  3,  1,  0, 44, 47,  0,...
    4,  0,  1,  1,  0, 50, 51,  0,...
    4,  0,  1,  1,  0, 52, 53,  0,...
    4,  1,  1,  1,  0, 54, 55,  0,...
    4,  1,  1,  1,  0, 56, 57,  0,...
    5,  0,  3,  1,  0, 44, 47,  0,...
    5,  0,  1,  1,  0, 50, 51,  0,...
    5,  0,  1,  1,  0, 52, 53,  0,...
    5,  1,  1,  1,  0, 54, 55,  0,...
    5,  1,  1,  1,  0, 56, 57,  0];
  env = [
    0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, ...
    0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, ...
    0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, ...
    0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, ...
    0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00, ...
   -2.93097845e+00, -2.16411436e-01,  0.00000000e+00,  0.00000000e+00, ...
   -3.65521976e+00,  1.44092184e+00,  0.00000000e+00,  0.00000000e+00, ...
   -1.13322529e+00,  7.69345300e-02,  0.00000000e+00,  0.00000000e+00, ...
    2.55231135e+00,  2.10645881e-01,  0.00000000e+00,  0.00000000e+00, ...
    3.17549200e+00, -7.06268132e-01, -1.43347254e+00,  0.00000000e+00, ...
    3.17549200e+00, -7.06268132e-01,  1.43347254e+00,  0.00000000e+00, ...
    1.30100000e+01,  1.96200000e+00,  4.44600000e-01,  5.79764064e-01, ...
    9.83419580e-01,  1.11930215e+00,  1.22000000e-01,  5.21536727e-01, ...
    2.97400000e-02,  1.80934235e-01,  7.27000000e-01,  1.95840453e+00, ...
    1.41000000e-01,  2.52062526e-01,  1.17200000e+04,  1.75900000e+03, ...
    4.00800000e+02,  1.13700000e+02,  3.70300000e+01,  1.32700000e+01, ...
    5.02500000e+00,  1.01300000e+00,  2.01954000e+00,  3.75176051e+00, ...
    6.29675240e+00,  9.21467523e+00,  1.07298897e+01,  7.87818279e+00, ...
    2.29637501e+00,  3.94147511e-02, -8.94856161e-01, -1.70329657e+00, ...
   -2.78735973e+00, -4.44591655e+00, -5.28622990e+00, -5.71024998e+00, ...
   -1.94898449e+00,  2.79438767e+00,  3.02300000e-01,  1.03001520e+00, ...
    7.89600000e-02,  3.76331333e-01,  1.77000000e+01,  3.85400000e+00, ...
    1.04600000e+00,  6.63856461e+00,  5.25433167e+00,  2.28748174e+00, ...
    2.75300000e-01,  5.81757722e-01,  6.85600000e-02,  1.02346478e-01, ...
    1.18500000e+00,  3.51185438e+00,  3.32000000e-01,  3.78896883e-01];
  natm = numel(atm)/6;
  nbas = numel(bas)/8;
  nao = ao_loc(shls_slice(2)+1) - ao_loc(shls_slice(1)+1);
  triuidx = 1:Norb^2;
  triflag = true(Norb,Norb);
  triflag = triu(triflag);
  triuidx = triuidx(triflag(:));
end

if 1 % 2025 version
  ngrids = n1*n2*n3;
  xyz = [x(:), y(:), z(:)]';
  xyz = xyz(:);
  fi = zeros([Norb ngrids]); 
  fi = GTOval_h2o_dimer_aug_ccpvdz_mwrap_mex(ngrids, shls_slice, ao_loc, fi, xyz, non0tab, atm, natm, bas, nbas, env);
  fi = reshape(fi,[Norb ngrids]);
  
  if 0 % naive version
    ftmp = zeros([ngrids nd]); % need to be changed later
    for k=1:ngrids
      fij = fi(:,k).*fi(:,k)';
      ftmp(k,:) = fij(triuidx);
    end
    f = reshape(ftmp,[n1 n2 n3 nd]);
  end

  % these will dominate due to ^2 scaling...
  firs = reshape(fi, [Norb, 1, ngrids]);
  fij_full = pagemtimes(firs, permute(firs, [2, 1, 3])); 
  ftmp = reshape(fij_full,[Norb^2 ngrids]);
  ftmp = ftmp(triuidx,:)';
  f = reshape(ftmp,[n1 n2 n3 nd]);
end  
