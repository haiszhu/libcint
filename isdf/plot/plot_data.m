%
%
%

molname = 'h2o_dimer';
% molname = 'nh3_dimer';

basmod = 'cc-pvdz.dat';
% basmod = 'aug-cc-pvdz.dat';
% basmod = 'cc-pvtz.dat';
% basmod = 'aug-cc-pvtz.dat';

[~, basname, ~] = fileparts(basmod);
mat_filename = [molname '_' basname];
load(mat_filename)

figure(1), clf
loglog(epsvals, relerrs, 'o-', 'LineWidth', 2)
hold on
loglog(epsvals, relerr2s, '^-', 'LineWidth', 2)
loglog(epsvals, relerr2ups, 's-', 'LineWidth', 2)

legend({ ...
    '$\max_{\tilde{T}} E_{\infty}^{(\phi_i)}$ on order $k$ grid', ...
    '$\max_{\tilde{T}} E_{\infty}^{(\rho_{ij})}$ on order $k$ grid', ...
    '$\max_{T} E_{\infty}^{(\rho_{ij})}$ on order $1.5\times k$ grid'}, ...
    'Interpreter','latex','FontSize',14,'Location','best')

xlabel('$\varepsilon$ using $L^2$ norm', 'Interpreter','latex','FontSize',14)
ylabel('achieved $E_{\infty}$', 'Interpreter','latex','FontSize',14)

if strcmp(molname, 'nh3_dimer') && strcmp(basmod, 'aug-cc-pvtz.dat')
title('NH$_3$ dimer, aug-cc-pvtz', 'Interpreter','latex','FontSize',14)
end
if strcmp(molname, 'nh3_dimer') && strcmp(basmod, 'cc-pvtz.dat')
title('NH$_3$ dimer, cc-pvtz', 'Interpreter','latex','FontSize',14)
end
if strcmp(molname, 'nh3_dimer') && strcmp(basmod, 'aug-cc-pvdz.dat')
title('NH$_3$ dimer, aug-cc-pvdz', 'Interpreter','latex','FontSize',14)
end
if strcmp(molname, 'nh3_dimer') && strcmp(basmod, 'cc-pvdz.dat')
title('NH$_3$ dimer, cc-pvdz', 'Interpreter','latex','FontSize',14)
end

if strcmp(molname, 'h2o_dimer') && strcmp(basmod, 'aug-cc-pvtz.dat')
title('H$_2$O dimer, aug-cc-pvtz', 'Interpreter','latex','FontSize',14)
end
if strcmp(molname, 'h2o_dimer') && strcmp(basmod, 'cc-pvtz.dat')
title('H$_2$O dimer, cc-pvtz', 'Interpreter','latex','FontSize',14)
end
if strcmp(molname, 'h2o_dimer') && strcmp(basmod, 'aug-cc-pvdz.dat')
title('H$_2$O dimer, aug-cc-pvdz', 'Interpreter','latex','FontSize',14)
end
if strcmp(molname, 'h2o_dimer') && strcmp(basmod, 'cc-pvdz.dat')
title('H$_2$O dimer, cc-pvdz', 'Interpreter','latex','FontSize',14)
end

grid on

keyboard