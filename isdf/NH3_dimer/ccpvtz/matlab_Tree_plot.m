% 
%
%

addpath('../../../')
addpath('../../../utils/')
addpath('../../../treefun/')

clear all

% setup
% isdf_base_path = '/mnt/home/cyeh/ceph/papers/isdf_adaptive/H2O_dimer/ccpvdz/isdf_adap/';
bdmk_exec = '../../../utils/f/int2-bdmk-mlscf';
treefun_order = 7;
treefun_eps = 1e-04; 
isdf_eps = 1e-3;
nd = 290;

% treefun_order = 4;
% treefun_eps = 1e-03; 
% nd = 80;

%%% resolve tree on cgto^2
rad = 7;
geom = sprintf([ ...
    'N  -1.578718  -0.046611   0.000000.\n',...
    'H  -2.158621   0.136396  -0.809565.\n',...
    'H  -2.158621   0.136396   0.809565.\n',...
    'H  -0.849471   0.658193   0.000000.\n',...
    'N   1.578718   0.046611   0.000000.\n',...
    'H   2.158621  -0.136396  -0.809565.\n',...
    'H   0.849471  -0.658193   0.000000.\n',...
    'H   2.158621  -0.136396   0.809565\n']),
molname = 'nh3_dimer';
basmod = 'cc-pvtz.dat';
basis = fullfile(fileparts(mfilename('fullpath')), '../../../basis', basmod);
mol = gto(geom,basis);
eval_name = 'GTOval_sph';
opts = struct('balance',true,...
              'tol',treefun_eps, ...
              'checkpts',mol.checkpts, ... 
              'ifcoeffs',false);
func = @(x,y,z) mol.eval_gto(eval_name, cat(4,x,y,z));
% func = @(x,y,z) mol.eval_gto2(eval_name, cat(4,x,y,z));
disp("=========Start treefun=======");
tic
f = treefun3(func,[-rad rad -rad rad -rad rad],treefun_order,opts);
time = toc;
disp("    treefun time is : " + time + " seconds");
disp("    treefun order is : " + treefun_order);
disp("    treefun num of boxes is : " + size(f.domain,2));
disp("=========End treefun=======");
disp("    ");
% figure(1),clf,plot(f,func);

% keyboard

nplotpts = 51;
[xx, yy, zz] = meshgrid(linspace(f.domain(1), f.domain(2), nplotpts), ...
                        linspace(f.domain(3), f.domain(4), nplotpts), ...
                        linspace(f.domain(5), f.domain(6), nplotpts));
tmpval = func(0,0,0);
nd = numel(tmpval);
v = zeros(size(xx));
tmpvals = func(xx,yy,zz);
for k = 1:nd
  v = v + tmpvals(:,:,:,k);
end
row = nplotpts;
col = nplotpts;
tube = 1;
xslice = xx(row,col,tube); 
yslice = yy(row,col,tube); 
zslice = zz(row,col,tube); 
yslice = [];
%
myLineWidth = 0.25;
%
figure(1),clf,
% slice(xx, yy, zz, v, xslice, yslice, zslice)
% shading interp
% % colorbar
% hold on
ids = leaves(f);
% figure(1),clf,
% xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids) ; nan(1, length(ids)); ... % bottom & top
%          f.domain([2 2 2 2 1 1], ids) ; nan(1, length(ids));]; % the rest to complete the box
% ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids) ; nan(1, length(ids)); ...
%          f.domain([3 3 4 4 4 4], ids) ; nan(1, length(ids));];
% zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids) ; nan(1, length(ids)); ...
%          f.domain([5 6 6 5 5 6], ids) ; nan(1, length(ids));];
% line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', 1, 'Color', [0 0 0]+0.85)
% axis equal, hold on
%
height = f.height(1);
ids1 = [];
for j=1:numel(ids)
  cenj = 1/2*(f.domain(1:2:5,ids(j))+f.domain(2:2:6,ids(j)));
  if (f.level(ids(j)) < height-4) && (-cenj(2)+cenj(3)<1 || cenj(1)>5  )
  % if (f.level(ids(j)) < height-3)
    ids1 = [ids1 ids(j)];
  end
end
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids1) ; nan(1, length(ids1)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids1) ; nan(1, length(ids1));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids1) ; nan(1, length(ids1)); ...
         f.domain([3 3 4 4 4 4], ids1) ; nan(1, length(ids1));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids1) ; nan(1, length(ids1)); ...
         f.domain([5 6 6 5 5 6], ids1) ; nan(1, length(ids1));];
% line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', 1, 'Color', [0 0 0]+0.85)
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth, 'Color', 'k')
axis equal, hold on

ids1color = [0.85 0.85 0.85];
for j=1:numel(ids1)
  dom = f.domain(:, ids1(j));
  drawboxfaces(dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), ids1color, 0.5);
end

ids2 = [];
for j=1:numel(ids)
  cenj = 1/2*(f.domain(1:2:5,ids(j))+f.domain(2:2:6,ids(j)));
  if (f.level(ids(j)) == height-4) && (-cenj(2)+cenj(3)<0.1 || cenj(1)>5  )
  % if (f.level(ids(j)) == height-3)
    ids2 = [ids2 ids(j)];
  end
end
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids2) ; nan(1, length(ids2)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids2) ; nan(1, length(ids2));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids2) ; nan(1, length(ids2)); ...
         f.domain([3 3 4 4 4 4], ids2) ; nan(1, length(ids2));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids2) ; nan(1, length(ids2)); ...
         f.domain([5 6 6 5 5 6], ids2) ; nan(1, length(ids2));];
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth, 'Color', 'm')
axis equal, hold on

ids2color = [1.0 0.7 1.0];
for j=1:numel(ids2)
  dom = f.domain(:, ids2(j));
  drawboxfaces(dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), ids2color, 0.5);
end

ids3 = [];
for j=1:numel(ids)
  cenj = 1/2*(f.domain(1:2:5,ids(j))+f.domain(2:2:6,ids(j)));
  if (f.level(ids(j)) == height-3) && -cenj(2)+cenj(3)<0.25
  % if (f.level(ids(j)) == height-2)
    ids3 = [ids3 ids(j)];
  end
end
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids3) ; nan(1, length(ids3)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids3) ; nan(1, length(ids3));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids3) ; nan(1, length(ids3)); ...
         f.domain([3 3 4 4 4 4], ids3) ; nan(1, length(ids3));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids3) ; nan(1, length(ids3)); ...
         f.domain([5 6 6 5 5 6], ids3) ; nan(1, length(ids3));];
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth, 'Color', 'b')
axis equal, hold on

ids3color = [0.6 0.8 1.0];
for j=1:numel(ids3)
  dom = f.domain(:, ids3(j));
  drawboxfaces(dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), ids3color, 0.5);
end

ids4 = [];
for j=1:numel(ids)
  cenj = 1/2*(f.domain(1:2:5,ids(j))+f.domain(2:2:6,ids(j)));
  if (f.level(ids(j)) == height-2) && (-cenj(2)+cenj(3)<0.1 || (cenj(2)>-0.2 && cenj(3)>1))
  % if (f.level(ids(j)) == height-1)
    ids4 = [ids4 ids(j)];
  end
end
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids4) ; nan(1, length(ids4)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids4) ; nan(1, length(ids4));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids4) ; nan(1, length(ids4)); ...
         f.domain([3 3 4 4 4 4], ids4) ; nan(1, length(ids4));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids4) ; nan(1, length(ids4)); ...
         f.domain([5 6 6 5 5 6], ids4) ; nan(1, length(ids4));];
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth/2, 'Color', 'g')
axis equal, hold on

ids4color = [0.7 1.0 0.7];
for j=1:numel(ids4)
  dom = f.domain(:, ids4(j));
  drawboxfaces(dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), ids4color, 0.5);
end

ids5 = [];
for j=1:numel(ids)
  cenj = 1/2*(f.domain(1:2:5,ids(j))+f.domain(2:2:6,ids(j)));
  % if (f.level(ids(j)) == height) && -cenj(1)-cenj(2)+cenj(3)<0
  if (f.level(ids(j)) == height-1) && (-cenj(2)+cenj(3)<0.025)
    ids5 = [ids5 ids(j)];
  end
end
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids5) ; nan(1, length(ids5)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids5) ; nan(1, length(ids5));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids5) ; nan(1, length(ids5)); ...
         f.domain([3 3 4 4 4 4], ids5) ; nan(1, length(ids5));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids5) ; nan(1, length(ids5)); ...
         f.domain([5 6 6 5 5 6], ids5) ; nan(1, length(ids5));];
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth/2, 'Color', 'r')
axis equal tight, hold on

ids5color = [1.0 0.6 0.6];
for j=1:numel(ids5)
  dom = f.domain(:, ids5(j));
  drawboxfaces(dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), ids5color, 0.5);
end

ids6 = [];
for j=1:numel(ids)
  cenj = 1/2*(f.domain(1:2:5,ids(j))+f.domain(2:2:6,ids(j)));
  % if (f.level(ids(j)) == height) && -cenj(1)-cenj(2)+cenj(3)<0
  if (f.level(ids(j)) == height)
    ids6 = [ids6 ids(j)];
  end
end
xdata = [f.domain([1 2 2 1 1 1 2 2 1 1], ids6) ; nan(1, length(ids6)); ... % bottom & top
         f.domain([2 2 2 2 1 1], ids6) ; nan(1, length(ids6));]; % the rest to complete the box
ydata = [f.domain([3 3 4 4 3 3 3 4 4 3], ids6) ; nan(1, length(ids6)); ...
         f.domain([3 3 4 4 4 4], ids6) ; nan(1, length(ids6));];
zdata = [f.domain([5 5 5 5 5 6 6 6 6 6], ids6) ; nan(1, length(ids6)); ...
         f.domain([5 6 6 5 5 6], ids6) ; nan(1, length(ids6));];
% line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth, 'Color', 'c')
line('XData', xdata(:), 'YData', ydata(:), 'ZData', zdata(:), 'LineWidth', myLineWidth/2, 'Color', [0.25 0.25 0.25])
axis equal tight, hold on

ids6color = [0.6 0.6 0.6];
for j=1:numel(ids6)
  dom = f.domain(:, ids6(j));
  drawboxfaces(dom(1), dom(2), dom(3), dom(4), dom(5), dom(6), ids6color, 0.5);
end

axis off
view(-45,45)
% title('Octree for NH$_3$ dimer', 'Interpreter','latex','FontSize',14)
% title('Octree for (NH_3)_2', 'FontSize',14)

keyboard

exportgraphics(gcf, 'octree.pdf','contenttype','vector');
% view(0,0)


function drawboxfaces(xmin, xmax, ymin, ymax, zmin, zmax, face_color, face_alpha)
v = [xmin ymin zmin;
     xmax ymin zmin;
     xmax ymax zmin;
     xmin ymax zmin;
     xmin ymin zmax;
     xmax ymin zmax;
     xmax ymax zmax;
     xmin ymax zmax];
faces = [1 2 3 4;  
         5 6 7 8;
         1 2 6 5;
         2 3 7 6;
         3 4 8 7;
         4 1 5 8];
for f = 1:6
  patch('Vertices', v, 'Faces', faces(f,:), 'FaceColor', face_color, 'EdgeColor', 'none', 'FaceAlpha', face_alpha);
end
end
