function [relerr,maxval] = treedata_resample_l2error(f,func,p)

%
ids = leaves(f);
[~,~,~,nd] = size(func(0,0,0));
x0 = (1-cos(pi*(2*(1:p)'-1)/(2*p)))/2;
[xx0, yy0, zz0] = ndgrid(x0);
% slightly different to verify results
nrefpts = p+1;
x0 = (1-cos(pi*(2*(1:nrefpts)'-1)/(2*nrefpts)))/2;
[xxx0, yyy0, zzz0] = ndgrid(x0);
x = 2*(x0-1/2);
l = floor(nrefpts/2)+1;
v = [2*exp(1i*pi*(0:nrefpts-l)/nrefpts)./(1-4*(0:nrefpts-l).^2)  zeros(1,l)];
w0 = real(ifft(v(1:nrefpts) + conj(v(nrefpts+1:-1:2))))';
[wx0, wy0, wz0] = ndgrid(w0/2); 
ww0 = wx0.*wy0.*wz0;
Eval = ones(nrefpts, p);
Eval(:,2) = x;
for k=3:p
    Eval(:,k) = 2*x.*Eval(:,k-1)-Eval(:,k-2);
end
errvalsl2s = zeros(nd,1);
maxvalsl2s = zeros(nd,1);
for k=1:numel(ids)
  % k-th leaf box
  idk = ids(k);
  domaink = f.domain(:,idk);
  % chebyshev points
  sclx = diff(domaink(1:2));
  scly = diff(domaink(3:4));
  sclz = diff(domaink(5:6));
  xx = sclx*xx0 + domaink(1); 
  yy = scly*yy0 + domaink(3); 
  zz = sclz*zz0 + domaink(5); 
  %
  F = func(xx,yy,zz);
  coeffsk = treefun3.vals2coeffs(F);

  % uniform points
  xxx = sclx*xxx0 + domaink(1); 
  yyy = scly*yyy0 + domaink(3); 
  zzz = sclz*zzz0 + domaink(5);
  %===================== order p grid =========================
  F = func(xxx,yyy,zzz);
  tmp1 = permute(tensorprod(Eval,coeffsk,2,1),[2 3 1 4]);
  tmp2 = permute(tensorprod(Eval,tmp1,2,1),[2 3 1 4]);
  vals = squeeze(permute(tensorprod(Eval,tmp2,2,1),[2 3 1 4]));

  %
  errvalsl2s = errvalsl2s + squeeze(sum((sclx*scly*sclz)*(F-vals).^2.*ww0,[1 2 3]));
  maxvalsl2s = maxvalsl2s + squeeze(sum((sclx*scly*sclz)*(F).^2.*ww0,[1 2 3]));
end
maxval_nd = sqrt(maxvalsl2s);
relerr_nd = sqrt(errvalsl2s)./maxval_nd;
relerr = max(relerr_nd);
maxval = max(maxval_nd);

end