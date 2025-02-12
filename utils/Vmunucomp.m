function Vmunu = Vmunucomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
                            ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                            nboxes,nlevels,ltree,itree,iptr,centers,boxsize)

%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
nd = size(fvals,2);
phi_kl = zeros(nd,npbox,nboxes); 
for k = 1:nd
  phi_kl(k,:,:) = reshape(fvals(:,k),[npbox nboxes]);
end
pot = zeros(nd,npbox,nboxes);
% pot = bdmk_wrap(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
%                  nboxes,nlevels,ltree,itree,iptr,centers,boxsize,phi_kl, ...
%                  pot);
pot = bdmk_wrap_mex(nd,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                 nboxes,nlevels,ltree,itree,iptr,centers,boxsize,phi_kl, ...
                 pot);
if 0 % wrap up
  nhess = ndim*(ndim+1)/2;
  nd = size(fvals,2);
  phi_kl = zeros(nd,npbox,nboxes); 
  for k = 1:nd
    phi_kl(k,:,:) = reshape(fvals(:,k),[npbox nboxes]);
  end
  pot = zeros(nd,npbox,nboxes);
  nd0 = 10;
  nchnk = floor(nd/nd0);
  disp("=========Start BDMK=======");
  disp("nchnk is : " + nchnk );
  % for jchnk = 1:1
  for jchnk = 1:nchnk
    jidxv = (jchnk-1)*nd0 + (1:nd0);
    phi_kl0 = phi_kl(jidxv,:,:);
    pot0=zeros(nd0,npbox,nboxes);
    grad0=zeros(nd0,ndim,npbox,nboxes);
    hess0=zeros(nd0,nhess,npbox,nboxes);
    tic
    ifpgh=1;
    ifpghtarg=0;
    ntarg = 100;
    targs=zeros(ndim,ntarg); pote=zeros(nd0,ntarg);
    grade=zeros(nd0,ndim,ntarg); hesse=zeros(nd0,nhess,ntarg);
    timeinfo = zeros(20,1);
    [pot0,grad0,hess0,pote,grade,hesse] = ...
      bdmk_mex(nd0,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
             nboxes,nlevels,ltree,itree,iptr,centers,boxsize,phi_kl0, ...
             ifpgh,pot0,grad0,hess0,ntarg,targs, ...
             ifpghtarg,pote,grade,hesse,timeinfo);
    pot(jidxv,:,:) = pot0;
    disp("jchnk is : " + jchnk + " ( total nchnk is : " + nchnk + " )");
  end
  jchnk = nchnk+1;
  jidxv = ((jchnk-1)*nd0+1):nd;
  nd0 = numel(jidxv);
  phi_kl0 = phi_kl(jidxv,:,:);
  pot0=zeros(nd0,npbox,nboxes);
  grad0=zeros(nd0,ndim,npbox,nboxes);
  hess0=zeros(nd0,nhess,npbox,nboxes);
  ifpgh=1;
  ifpghtarg=0;
  ntarg = 100;
  targs=zeros(ndim,ntarg); pote=zeros(nd0,ntarg);
  grade=zeros(nd0,ndim,ntarg); hesse=zeros(nd0,nhess,ntarg);
  timeinfo = zeros(20,1);
  [pot0,grad0,hess0,pote,grade,hesse] = ...
  bdmk_mex(nd0,ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
         nboxes,nlevels,ltree,itree,iptr,centers,boxsize,phi_kl0, ...
         ifpgh,pot0,grad0,hess0,ntarg,targs, ...
         ifpghtarg,pote,grade,hesse,timeinfo);
  pot(jidxv,:,:) = pot0;
  disp("=========End BDMK=======");
end

%
phi_ij_leaf = zeros(nd,npbox,nleafbox); % leaf box information
potleaf = zeros(nd,npbox,nleafbox);
jbox = 0;
for ilev = 0:nlevels
  bs = boxsize(ilev+1); % Adjusted for MATLAB indexing (1-based)
  ifirstbox = itree(2 * ilev + 1);
  ilastbox = itree(2 * ilev + 2);
  nbloc = ilastbox - ifirstbox + 1;
  for i = 1:nbloc
    ibox = ifirstbox + i - 1;
    if itree(iptr(4) + ibox - 1) == 0
      jbox = jbox + 1;
      for ell = 1:npbox
        for j = 1:nd
          phi_ij_leaf(j, ell, jbox) = phi_kl(j, ell, ibox);
          potleaf(j, ell, jbox) = pot(j, ell, ibox);
        end
      end
    end
  end
end

%% compute smooth volume integral
disp("=========Start Smooth Integral=======");
% nd = nd0;
potleaf = reshape(potleaf,nd,[])/ratio^2; 
phi_ij_leaf = reshape(phi_ij_leaf,nd,[]);
wtsleaf = wtsleaf(:)';
phi_ij_leaf = phi_ij_leaf.*wtsleaf;
disp("nd is : " + nd );
Vmunu = phi_ij_leaf*potleaf';
Vmunu = Vmunu/ratio^3;
if 0 % old version pre 2025
  Vmunu = zeros(nd,nd);
  % ell = 1;
  for ell = 1:nd
    % (2*L)^2 * potleaf = (2*L)^2 * \int_{-1/2}^{1/2} (phi_k*phi_l)/|r' - r| dV
    potleaf_ell = potleaf(ell,:,:)/ratio^2; 
    for j = 1:nd
      % (2*L)^3 * \int_{-1/2}^{1/2} phi_{ij} * pot_{kl} dV
      Vmunu(j,ell) = sum(squeeze(phi_ij_leaf(j,:,:).*potleaf_ell).*wtsleaf(:,:),'all')/ratio^3;
    end
  end
end
disp("=========End Smooth Integral=======");

end