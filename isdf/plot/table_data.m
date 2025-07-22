% table_data.m
% summarize all basis results for one molecule

clear; clc;

% molecule to process
molname = 'h2o_dimer';
% molname = 'nh3_dimer';

% list of basis sets to process
% basis_list = {'cc-pvdz.dat', 'aug-cc-pvdz.dat', 'cc-pvtz.dat', 'aug-cc-pvtz.dat'};
basis_list = {'aug-cc-pvtz.dat'};

% initialize output table
T = table();

for ibasis = 1:length(basis_list)
    basmod = basis_list{ibasis};
    [~, basname, ~] = fileparts(basmod);
    mat_filename = [molname '_' basname '.mat'];
    
    % load data
    fprintf('Loading %s...\n', mat_filename);
    load(mat_filename);
    
    % assume epsvals are consistent across files
    if ibasis == 1
        epsvals_all = epsvals(:); % make sure column
    end
    
    % add columns to table
    T.(['relerrs_' basname]) = relerrs(:);
    T.(['relerr2s_' basname]) = relerr2s(:);
    T.(['relerr2ups_' basname]) = relerr2ups(:);
end

% add epsvals as first column
T = addvars(T, epsvals_all, 'Before', 1, 'NewVariableNames', 'epsvals');

% display table
disp(T)

keyboard
% optional: save as csv
output_csv = [molname '_summary.csv'];
writetable(T, output_csv);
fprintf('Table saved to %s\n', output_csv);
