function relerr = treedata_error(f,func)
%
%

ids = leaves(f);
[p,~,~,nd] = size(f.coeffs{ids(1)});
nrefpts = p;
[xxx0, yyy0, zzz0] = ndgrid(linspace(0, 1, nrefpts));
x = linspace(-1, 1, nrefpts).';
Eval = ones(nrefpts, p);
Eval(:,2) = x;
for k=3:p
    Eval(:,k) = 2*x.*Eval(:,k-1)-Eval(:,k-2);
end
err = zeros(numel(ids),1);
maxval = zeros(numel(ids),1);
for k=1:numel(ids)
  % k-th leaf box
  idk = ids(k);
  domaink = f.domain(:,idk);
  levelk = f.level(idk);
  heightk = f.height(idk);
  parentk = f.parent(idk);
  childrenk = f.children(:,idk); % all 0
  coeffsk = f.coeffs{idk};
  % uniform points
  sclx = diff(domaink(1:2));
  scly = diff(domaink(3:4));
  sclz = diff(domaink(5:6));
  h = sclx;
  xxx = sclx*xxx0 + domaink(1); 
  yyy = scly*yyy0 + domaink(3); 
  zzz = sclz*zzz0 + domaink(5);
  %===================== order p grid =========================
  F = func(xxx,yyy,zzz);
  tmp1 = permute(tensorprod(Eval,coeffsk,2,1),[2 3 1 4]);
  tmp2 = permute(tensorprod(Eval,tmp1,2,1),[2 3 1 4]);
  vals = squeeze(permute(tensorprod(Eval,tmp2,2,1),[2 3 1 4]));
  err(k) = max(abs(F(:)-vals(:)));
  maxval(k) = max(abs(F(:)));

end
relerr = max(err)/max(maxval);

end