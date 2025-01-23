function [ao] = GTOval_uracil_dimer_ccpvdz_mwrap_mex( ngrids, shls_slice, ao_loc, ao, coord, non0table, atm, natm, bas, nbas, env)
naoloc = shls_slice(2) - shls_slice(1) + 1;
nao = ao_loc(shls_slice(2)+1) - ao_loc(shls_slice(1)+1);
ntmpao = ngrids*nao;
ndim = 3;
ndimngrids = numel(coord);
nnon0tab = numel(non0table);
ATM_SLOTS = 6;
BAS_SLOTS = 8;
ntmpatm = natm*ATM_SLOTS;
ntmpbas = nbas*BAS_SLOTS;
nenv = numel(env);
mex_id_ = 'GTOval_uracil_dimer_ccpvdz_mwrap(i int[x], i int[x], i int[x], io double[x], i double[x], i int[x], i int[x], i int[x], i int[x], i int[x], i double[x])';
[ao] = gateway(mex_id_, ngrids, shls_slice, ao_loc, ao, coord, non0table, atm, natm, bas, nbas, env, 1, 2, naoloc, ntmpao, ndimngrids, nnon0tab, ntmpatm, 1, ntmpbas, 1, nenv);
end
