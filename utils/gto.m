function mol = gto(atom,basis)
% some functionality of pyscf gto
% 
%
% 02/10/25, Hai

% constants in Python code
CHARGE_OF   = 1; % indexing +1
PTR_COORD   = 2;
NUC_MOD_OF  = 3;
PTR_ZETA    = 4;
ATM_SLOTS   = 6; 
ATOM_OF     = 1; % indexing +1
PTR_ENV_START   = 20;
%
NUC_POINT   = 1;

%
BOHR2ANG = 0.529177210903;

%
unit = 'Angstrom';

%
geom = atom;

% split and store
lines = splitlines(geom);
atoms_list = {};
for i = 1:numel(lines)
  line = strtrim(lines{i});
  if isempty(line)
    continue
  end
  parts = strsplit(line);
  if numel(parts) < 4
    continue  % 
  end
  symb = parts{1};
  x = str2double(parts{2});
  y = str2double(parts{3});
  z = str2double(parts{4});
  % 
  atoms_list{i,1} = symb;
  atoms_list{i,2} = [x, y, z];
end

% format atom
nAtoms = size(atoms_list, 1);
scale = 1.0;
if strcmpi(unit,'Angstrom')
  scale = 1.0 / BOHR2ANG;
end
for i = 1:nAtoms
  symb = atoms_list{i,1};
  coords = atoms_list{i,2}; % in Angstrom
  coords_bohr = coords * scale;
  atoms_list_bohr{i,1} = symb;
  atoms_list_bohr{i,2} = coords_bohr;  % now in Bohr
end

% load(basisfile, symb, optimize...)
all_symbols = atoms_list_bohr(:,1);    
unique_elems = unique(all_symbols);   

% create an empty containers.Map for the final basis_dict
basis_dict = containers.Map('KeyType','char', 'ValueType','any');
for ielem = 1:numel(unique_elems)
  symb = unique_elems{ielem};
  % load & store
  basis_dict(symb) = loadNWChemBasis(basis, symb);
end

%
pre_env = zeros(1,PTR_ENV_START);
% ELEMENTS_PROTON = containers.Map({'N','H'}, [7, 1]);
ELEMENTS_PROTON = buildElementsProtonMap();

%
atm = zeros(nAtoms, ATM_SLOTS);
env = [pre_env];
ptr_env = numel(pre_env);
for ia = 1:nAtoms
  % 
  symb   = atoms_list_bohr{ia, 1};
  coords = atoms_list_bohr{ia, 2}; 

  % 
  nuc_charge = ELEMENTS_PROTON(symb);  

  % 
  zeta = 0.0;
  nuclear_model = NUC_POINT;

  % 
  env0 = [coords, zeta];
  len_env0 = numel(env0);  

  % 
  atm0 = zeros(1, ATM_SLOTS);
  atm0(CHARGE_OF)  = nuc_charge;
  atm0(PTR_COORD)  = ptr_env;         
  atm0(NUC_MOD_OF) = nuclear_model;   
  atm0(PTR_ZETA)   = ptr_env + 3;     

  % 
  atm(ia,:) = atm0;
  env = [env, env0];  
  ptr_env = ptr_env + len_env0; 
end

% initialize basis dictionary
if isa(basis_dict, 'struct')
  allSymbols = fieldnames(basis_dict);
  isMapInput = false;
else
  % if basis_dict is a containers.Map
  allSymbols = keys(basis_dict);
  isMapInput = true;
end
% 
basdic = containers.Map();
for iSym = 1:numel(allSymbols)
  symb = allSymbols{iSym};

  if isMapInput
    basis_add = basis_dict(symb);
  else
    basis_add = basis_dict.(symb);
  end

  ptr_env0 = ptr_env;
  atom_id = 0;
  bas0 = [];   
  env0 = [];
  for ib = 1:numel(basis_add)
    b = basis_add{ib};  
    angl = b{1};  
    kappa = 0;
    b_coefftmp = b(2:end);  
    b_coeff = cell2mat(b_coefftmp(:));  
    [~, idxSort] = sort(b_coeff(:,1), 'descend');
    b_coeff = b_coeff(idxSort,:);
    es = b_coeff(:,1);
    cs = b_coeff(:,2:end);  
    [nprim, nctr] = size(cs);

    n1 = double((angl*2 + 2) + 1) * 0.5;
    gaussian_intangles = gamma(n1) ./ (2 .* ((2.0 .* es).^n1));
    gto_normangles = 1.0 ./ sqrt(gaussian_intangles);
    for ip = 1:nprim
      cs(ip,:) = cs(ip,:) * gto_normangles(ip);
    end

    ee = es(:) + es(:).';
    n1 = ((angl*2 + 2) + 1) * 0.5;
    ee = gamma(n1) ./ (2 .* (ee.^n1));
    % double check
    scales = zeros(1, nctr);
    for ic = 1:nctr
      Cp = cs(:, ic);
      val = 0;
      for p = 1:nprim
        for q = 1:nprim
          val = val + Cp(p)* ee(p,q)* Cp(q);
        end
      end
      scales(ic) = 1.0 / sqrt(val);
    end
    for ic = 1:nctr
      cs(:,ic) = cs(:,ic) * scales(ic);
    end

    env0 = [env0, [ reshape(es,1,[]), reshape(cs.',1,[]) ]];  
    ptr_exp = ptr_env0;                
    ptr_coeff = ptr_exp + nprim;         
    ptr_env0 = ptr_coeff + nprim*nctr;  
    bas0 = [bas0; [atom_id, angl, nprim, nctr, kappa, ptr_exp, ptr_coeff, 0]]; 
  end

  env0 = env0;
  bas0 = bas0;

  ptr_env = ptr_env + numel(env0);
  basdic(symb) = bas0;
  env = [env, env0];
end

%
% [bas] = gatherBasForAtoms(basdic, atoms_list_bohr);
basBlocks = cell(nAtoms, 1);
for ia=1:nAtoms
  atomi = atoms_list_bohr(ia,:);
  symb = atomi{1};  % e.g. 'N' or 'H'
  b = basdic(symb);
  b(:, ATOM_OF) = ia-1;
  basBlocks{ia} = b;
end
bas = vertcat(basBlocks{:});

%
mol = [];
mol.atm = atm;
mol.bas = bas;
mol.env = env;
end
