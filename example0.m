% 
%
% 12/01/24 Hai


addpath('./utils/')
addpath('./treefun/')

clear all

% 
Order_all = [4 5 6 7 8 9 10];
Eps_all = [1e-03 1e-04 1e-05 1e-06 1e-07 1e-08 1e-09];
idfname{1} = 'example0_1e-3.h5'; % 0.0026
idfname{2} = 'example0_1e-4.h5';
idfname{3} = 'example0_1e-5.h5';
idfname{4} = 'example0_1e-6.h5';
idfname{5} = 'example0_1e-7.h5';
idfname{6} = 'example0_1e-8.h5';
idfname{7} = 'example0_1e-9.h5';

%
for iname = 1:numel(Eps_all)
  order = Order_all(iname);
  eps = Eps_all(iname);
  %
  %%% resolve tree on cgto^2
  func2 = @(x,y,z) cgto2func(x,y,z);
  opts = struct('balance',true,...
                'tol',eps, ...
                'checkpts',[0 0 0; 0 -0.757 0.757;0 0.587 0.587], ...
                'ifcoeffs',false);
  f = treefun3(func2,[-15 15 -15 15 -15 15],order,opts); 

  %%% treefun to bdmk
  Norb = 24; % 
  ndim = 3;
  ratio = 0.5/15; % from boxlen to 1
  ipoly = 0;
  [src,nleafbox,srcleaf,wtsleaf,norder,npbox,nboxes,nlevels,ltree,itree,iptr,centers,boxsize] = treefun2bdmk(f,ndim,ratio,ipoly);

  %%% eval cgto
  src0 = src/ratio;
  func = @(x,y,z) cgtofunc(x,y,z);
  fvals0 = squeeze(func(squeeze(src0(1,:,:)),squeeze(src0(2,:,:)),squeeze(src0(3,:,:))));
  fvals0 = permute(fvals0,[3 1 2]);
  fvals = fvals0;
  
  %%% compute V_ijkl
  nd = Norb*(Norb+1)/2;
  ikernel = 1;
  beta = 6.0d0;
  Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
                  ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                  nboxes,nlevels,ltree,itree,iptr,centers,boxsize);
  
  %%% save data
  save(idfname{iname}, 'Vijkl');
  
  % keyboard
end

keyboard