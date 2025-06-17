function [relerr,maxval] = treedata_resample_error(f,func,p)

%
ids = leaves(f);
[~,~,~,nd] = size(func(0,0,0));
x0 = (1-cos(pi*(2*(1:p)'-1)/(2*p)))/2;
[xx0, yy0, zz0] = ndgrid(x0);

nrefpts = p;
[xxx0, yyy0, zzz0] = ndgrid(linspace(0, 1, nrefpts));
x = linspace(-1, 1, nrefpts).';
Eval = ones(nrefpts, p);
Eval(:,2) = x;
for k=3:p
    Eval(:,k) = 2*x.*Eval(:,k-1)-Eval(:,k-2);
end
errvals = zeros(numel(ids),nd);
maxvals = zeros(numel(ids),nd);

%
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
  errvals(k,:) = max(reshape(abs(F-vals),[],nd),[],1);
  maxvals(k,:) = max(reshape(abs(F),[],nd),[],1);
end
maxval_nd = max(maxvals,[],1);
relerr_nd = max(errvals,[],1)./maxval_nd;
relerr = max(relerr_nd);
maxval = max(maxval_nd);

end