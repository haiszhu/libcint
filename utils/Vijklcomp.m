function Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
                            ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                            nboxes,nlevels,ltree,itree,iptr,centers,boxsize)
nhess = ndim*(ndim+1)/2;
%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
nd = Norb*(Norb+1)/2;
phi_kl = zeros(nd,npbox,nboxes); 
idx = 0;
for k = 1:Norb
  for ell = k:Norb
    idx = idx + 1;
    phi_kl(idx,:,:) = fvals(k,:,:).*fvals(ell,:,:);
  end
end
pot = zeros(nd,npbox,nboxes);
nd0 = 10;
nchnk = nd/nd0;
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
  time = toc;
end
% ntot=npbox*nleafbox+ntarg*ifpghtarg;
% X = [' ntotal= ',num2str(ntot)];
% disp(X)
% pps = ntot/time;
% X = [' speed in pps= ', num2str(pps)];
% disp(X)

xleaf = srcleaf(1,:,:); xleaf = xleaf(:);
yleaf = srcleaf(2,:,:); yleaf = yleaf(:);
zleaf = srcleaf(3,:,:); zleaf = zleaf(:);
phi_ij_leaf = zeros(nd,npbox,nleafbox); % leaf box information
potleaf = zeros(nd,npbox,nleafbox);
potnonleaf = zeros(nd,npbox,nboxes-nleafbox);
jbox = 0;
jbox2 = 0;
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
    else
      jbox2 = jbox2 + 1;
      for ell = 1:npbox
        for j = 1:nd
          potnonleaf(j, ell, jbox2) = pot(j, ell, ibox);
        end
      end
    end
  end
end

%% compute smooth volume integral
% nd = nd0;
Vijkl0 = zeros(nd,nd);
% ell = 1;
for ell = 1:nd
  % (2*L)^2 * potleaf = (2*L)^2 * \int_{-1/2}^{1/2} (phi_k*phi_l)/|r' - r| dV
  potleaf_ell = potleaf(ell,:,:)/ratio^2; 
  for j = 1:nd
    % (2*L)^3 * \int_{-1/2}^{1/2} phi_{ij} * pot_{kl} dV
    Vijkl0(j,ell) = sum(squeeze(phi_ij_leaf(j,:,:).*potleaf_ell).*wtsleaf(:,:),'all')/ratio^3;
  end
end

% Vijkl = Vijkl0;
% turn back into 4-tensor and save
Vijkltmp = zeros(Norb,Norb,nd);
for kl = 1:nd
  tmpidx = 0;
  for i = 1:Norb
    for j = i:Norb
      tmpidx = tmpidx + 1;
      Vijkltmp(i,j,kl) = Vijkl0(tmpidx,kl);
    end
  end
  for i = 1:Norb
    for j = 1:i-1
      Vijkltmp(i,j,kl) = Vijkltmp(j,i,kl);
    end
  end
end
Vijkl = zeros(Norb,Norb,Norb,Norb);
tmpidx = 0;
for k = 1:Norb
  for l = k:Norb
    tmpidx = tmpidx + 1;
    Vijkl(:,:,k,l) = Vijkltmp(:,:,tmpidx);
  end
end
for k = 1:Norb
  for l = 1:k-1
    Vijkl(:,:,k,l) = Vijkl(:,:,l,k);
  end
end
end
