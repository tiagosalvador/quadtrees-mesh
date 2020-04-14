%% Example script - Solves the Convex Envelope equation

%% Setup
max_depth = 6; % depth of the quadtree;
% it corresponds to a N by N uniform grid with N = 2^(depth-1)+1;
% Computational domain: [a,b]x[a,b]
a = -1;
b = 1;
% Function to determine the boundary
domain_function = @domain_conv_env;
% Exact Solution
uex = @uConvEnv;
% Obstacle function in Convex Envelope equation
g = @gConvEnv;

% Function to determine dtheta
dTh_choice = @(h) sqrt(h);

%% Solving Convex Envelope equation
[uM,uA,grid_parameters,mesh_parameters,gridR_parameters,meshR_parameters] = solveConvEnvAdaptivity(a,b,domain_function,max_depth,uex,g,dTh_choice);

%% Plotting solution
x = grid_parameters.x;
y = grid_parameters.y;
tri = delaunay(x,y);
figure;
trisurf(tri,x,y,uM)
title('Solution on standard mesh')
shading flat
errorM = norm(uex(x,y)-uM,Inf);

%% Plotting solution
xR = gridR_parameters.x;
yR = gridR_parameters.y;
tri = delaunay(xR,yR);
figure;
trisurf(tri,xR,yR,uA)
title('Solution on refined mesh')
shading flat
errorA = norm(uex(xR,yR)-uA,Inf);

%% Drawing mesh
figure;
plot(grid_parameters.x,grid_parameters.y,'.k',grid_parameters.xB,grid_parameters.yB,'.b')
title('Mesh')
axis([a b a b])

figure;
plot(gridR_parameters.x,gridR_parameters.y,'.k',gridR_parameters.xB,gridR_parameters.yB,'.b')
title('Refined Mesh')
axis([a b a b])
