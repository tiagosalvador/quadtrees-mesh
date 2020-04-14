%% Example script - Draws the quadtree meshes on different domains

%% Ellipse
max_depth = 5; % depth of the quadtree;
% it corresponds to a N by N uniform grid with N = 2^(depth-1)+1;
% Computational domain: [a,b]x[a,b]
a = -1;
b = 1;
% Function to determine the boundary
domain_function = @domain_circle;

% Function to determine dtheta
dTh_choice = @(h) 2*nthroot(h,3);
[grid_parameters,~,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice);
draw_mesh(a,b,grid_parameters,mesh_parameters,domain_function);
set(gca,'Visible','off')
% print(gcf, '-depsc2', 'circle.eps');

%% Ellipse
max_depth = 5; % depth of the quadtree;
% it corresponds to a N by N uniform grid with N = 2^(depth-1)+1;
% Computational domain: [a,b]x[a,b]
a = -1;
b = 1;
% Function to determine the boundary
domain_function = @domain_ellipse;

% Function to determine dtheta
dTh_choice = @(h) 2*nthroot(h,3);
[grid_parameters,~,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice);
draw_mesh(a,b,grid_parameters,mesh_parameters,domain_function);
set(gca,'Visible','off')
% print(gcf, '-depsc2', 'ellipse.eps');

%% Diammond Stretched
max_depth = 5; % depth of the quadtree;
% it corresponds to a N by N uniform grid with N = 2^(depth-1)+1;
% Computational domain: [a,b]x[a,b]
a = -1;
b = 1;
% Function to determine the boundary
domain_function = @domain_diamond_stretched;

% Function to determine dtheta
dTh_choice = @(h) 2*nthroot(h,3);
[grid_parameters,~,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice);
draw_mesh(a,b,grid_parameters,mesh_parameters,domain_function);
set(gca,'Visible','off')
% print(gcf, '-depsc2', 'dimamond_stretched.eps');

%% Clover
max_depth = 5; % depth of the quadtree;
% it corresponds to a N by N uniform grid with N = 2^(depth-1)+1;
% Computational domain: [a,b]x[a,b]
a = -1.5;
b = 1.5;
% Function to determine the boundary
domain_function = @domain_clover;

% Function to determine dtheta
dTh_choice = @(h) 2*nthroot(h,3);
[grid_parameters,~,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice);
draw_mesh(a,b,grid_parameters,mesh_parameters,domain_function);
set(gca,'Visible','off')
% print(gcf, '-depsc2', 'clover.eps');