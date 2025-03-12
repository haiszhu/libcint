function Vijkl = computeVijkl_mex(nd, Norb, idcoefs, Vmunu, Vijkl)
Norbtp1d2 = Norb*(Norb+1)/2;
Norb2 = Norb^2;
Vijkl = reshape(Vijkl,[Norb2 Norb2]);
mex_id_ = 'computeVijkl(i int[x], i int[x], i double[xx], i double[xx], io double[xx])';
[Vijkl] = helper_module(mex_id_, nd, Norb, idcoefs, Vmunu, Vijkl, 1, 1, nd, Norbtp1d2, nd, nd, Norb2, Norb2);
Vijkl = reshape(Vijkl,[Norb Norb Norb Norb]);
end