%% Example script - Solves the Monge-Ampere equation on an ellipse

%% Setup
max_depth = 5; % depth of the quadtree;
% it corresponds to a N by N uniform grid with N = 2^(depth-1)+1;
% Computational domain: [a,b]x[a,b]
a = -1;
b = 1;
% Function to determine the boundary
domain_function = @domain_ellipse; 
% Exact Solution
uex = @uMA;
% Right-hand side of MA equation
f = @fMA;

% Function to determine dtheta
dTh_choice = @(h) 2*nthroot(h,3);
% Function to determine the choice of epsilon for the filtered scheme
epsilon_choice = @(h,dTh) sqrt(h/dTh+dTh);

%% Solving Monge-Ampere equation
[uM,uF,grid_parameters,mesh_parameters] = solveMA(a,b,domain_function,max_depth,uex,f,dTh_choice,epsilon_choice);

%% Plotting solutions
x = grid_parameters.x;
y = grid_parameters.y;
tri = delaunay(x,y);
figure;
trisurf(tri,x,y,uM)
title('Solution Monotone Scheme')
shading flat
figure;
trisurf(tri,x,y,uF)
title('Solution Filtered Scheme')
shading flat

%% Compute errors
errorM = norm(uex(x,y)-uM,Inf);
errorF = norm(uex(x,y)-uF,Inf);

%% Drawing mesh
figure;
indicator_boundary = @(x,y) domain_function(x,y);
insidedomain = @(x,y) domain_function(x,y)<0;
draw_quadtree(grid_parameters,mesh_parameters,insidedomain);
[q1,q2] = meshgrid(linspace(-2,2,1000));
hold on
contour(q1,q2,indicator_boundary(q1,q2),[0 0],'LineColor','black','LineWidth',2)
x = grid_parameters.x;
y = grid_parameters.y;
xB = grid_parameters.xB;
yB = grid_parameters.yB;
plot(x,y,'ok',xB,yB,'ob')
title('Mesh')
axis([a b a b])
hold off