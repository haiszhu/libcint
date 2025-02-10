% matlab build mole
%
% 02/10/25 Hai

addpath('./utils/')

% geom
geom = sprintf([ ...
    'N  -1.578718  -0.046611   0.000000\n',...
    'H  -2.158621   0.136396  -0.809565\n',...
    'H  -2.158621   0.136396   0.809565\n',...
    'H  -0.849471   0.658193   0.000000\n',...
    'N   1.578718   0.046611   0.000000\n',...
    'H   2.158621  -0.136396  -0.809565\n',...
    'H   0.849471  -0.658193   0.000000\n',...
    'H   2.158621  -0.136396   0.809565\n']),

% mol
basmod = 'cc-pvdz.dat';
% basmod = 'aug-cc-pvdz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), 'basis', basmod);
mol = gto(geom,basis);

% atm, bas, env
mol.atm
mol.bas
reshape(mol.env,4,[])'
mol.ao_loc
mol.non0tab

keyboard

