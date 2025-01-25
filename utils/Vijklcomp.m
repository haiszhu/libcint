function Vijkl = Vijklcomp(Norb,ratio,fvals,nleafbox,srcleaf,wtsleaf,...
                            ndim,eps,ikernel,beta,ipoly,norder,npbox, ...
                            nboxes,nlevels,ltree,itree,iptr,centers,boxsize)
nhess = ndim*(ndim+1)/2;
%%% based on fvals for all cGTOs, compute phi_k(r')*phi_l(r'), and singular volume integral
nd = Norb*(Norb+1)/2;
phi_kl = zeros(nd,npbox*nboxes); 
idx = 0;
disp("=========Start phi_kl=======");
disp("nd is : " + nd );
fvals_rs = reshape(fvals,[Norb npbox*nboxes]);
for k = 1:Norb
  phi_kl(idx+(1:Norb-k+1),:) = fvals_rs(k,:).*fvals_rs(k:Norb,:);
  idx = idx + (Norb-k+1);
  disp("k in phi_kl is : " + k + " ( total Norb is : " + Norb + " )");
end
phi_kl = reshape(phi_kl,nd,npbox,nboxes);
if 0 % naive pre 2025
  phi_kl = zeros(nd,npbox,nboxes); 
  for k = 1:Norb
    for ell = k:Norb
      idx = idx + 1;
      phi_kl(idx,:,:) = fvals(k,:,:).*fvals(ell,:,:);
    end
  end
end
disp("=========End phi_kl=======");

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
  time = toc;
  disp("jchnk is : " + jchnk + " ( total nchnk is : " + nchnk + " )");
end
jchnk = nchnk + 1;
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
disp("=========Start Smooth Integral=======");
% nd = nd0;
potleaf = reshape(potleaf,nd,[])/ratio^2; 
phi_ij_leaf = reshape(phi_ij_leaf,nd,[]);
wtsleaf = wtsleaf(:)';
phi_ij_leaf = phi_ij_leaf.*wtsleaf;
disp("nd is : " + nd );
Vijkl0 = phi_ij_leaf*potleaf';
Vijkl0 = Vijkl0/ratio^3;
if 0 % old version pre 2025
  Vijkl0 = zeros(nd,nd);
  potleaf = reshape(potleaf,nd,[])/ratio^2; 
  phi_ij_leaf = reshape(phi_ij_leaf,nd,[]);
  wtsleaf = wtsleaf(:)';
  disp("nd is : " + nd );
  for ell = 1:nd
    % (2*L)^2 * potleaf = (2*L)^2 * \int_{-1/2}^{1/2} (phi_k*phi_l)/|r' - r| dV
    Vijkl0(:,ell) = phi_ij_leaf*(potleaf(ell,:).*wtsleaf)';
    disp("ell is : " + ell + " ( total nd is : " + nd + " )");
  end
  Vijkl0 = Vijkl0/ratio^3;
end
% Vijkl0 = zeros(nd,nd);
% % ell = 1;
% for ell = 1:nd
%   % (2*L)^2 * potleaf = (2*L)^2 * \int_{-1/2}^{1/2} (phi_k*phi_l)/|r' - r| dV
%   potleaf_ell = potleaf(ell,:,:)/ratio^2; 
%   for j = 1:nd
%     % (2*L)^3 * \int_{-1/2}^{1/2} phi_{ij} * pot_{kl} dV
%     Vijkl0(j,ell) = sum(squeeze(phi_ij_leaf(j,:,:).*potleaf_ell).*wtsleaf(:,:),'all')/ratio^3;
%   end
% end
disp("=========End Smooth Integral=======");

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
