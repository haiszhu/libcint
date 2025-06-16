function g = treedata_resample(f,func,p)
% take tree for f as given
% resample coefficients to order p for some function handle
% 
% 06/16/25 Hai

%
g = f;

%
ids = leaves(f);
x0 = (1-cos(pi*(2*(1:p)'-1)/(2*p)))/2;
[xx0, yy0, zz0] = ndgrid(x0);

%
for k=1:numel(ids)
  % k-th leaf box
  idk = ids(k);
  domaink = f.domain(:,idk);
  % uniform points
  sclx = diff(domaink(1:2));
  scly = diff(domaink(3:4));
  sclz = diff(domaink(5:6));
  xx = sclx*xx0 + domaink(1); 
  yy = scly*yy0 + domaink(3); 
  zz = sclz*zz0 + domaink(5); 
  %
  F = func(xx,yy,zz);
  g.coeffs{idk} = treefun3.vals2coeffs(F);

end

end