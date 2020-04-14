function [uM,uF,grid_parameters,mesh_parameters] = solveMA(a,b,domain_function,max_depth,uex,f,dTh_choice,epsilon_choice)
% Solves the Monge-Ampere equation on [a,b]x[a,b] on a uniform grid of
% size n by n.
% Args:
%    domain_function: variable to determine the domain
%    max_depth: maximum depth of the initial quadtree
%    uex: exact solution to impose boundary condtions
%    f: right hand side of equation
%    dTh_choice: function to determine dTh
%    epsilon_choice: function to determine epsilon for the filtered scheme
% Returns:
%    uM: approximate solution from the monotone scheme
%    uF: approximate solution from the filtered scheme
%    grid_paramters: Matlab structure containing grid information
%    mesh_paramters: Matlab structure containing mesh information



%% Setup

indicator_boundary = @(x,y) domain_function(x,y);
insidedomain = @(x,y) domain_function(x,y)<0;

[grid_parameters,dTh,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice);
N = grid_parameters.N;
h = grid_parameters.h;
xgrid = grid_parameters.xgrid;
ygrid = grid_parameters.ygrid;
xgrid2x = grid_parameters.xgrid2x;
iBdy = grid_parameters.iBdy;
iInt = grid_parameters.iInt;
x = grid_parameters.x;
y = grid_parameters.y;
numTotGrid = grid_parameters.numTotGrid;
numInt = grid_parameters.numInt;
numBdy = grid_parameters.numBdy;
numTot = grid_parameters.numTot;

quadtree = mesh_parameters.quadtree;
squares_boundary = mesh_parameters.squares_boundary;
vertices_square = mesh_parameters.vertices_square;
squares_vertex = mesh_parameters.squares_vertex;


%% Finding neighbours for monotone scheme
r = h*(1+cos(dTh/2)*cot(dTh/2)+sin(dTh/2));
numFr = round(0.5*pi/sqrt(dTh));
numDir = 2*numFr;
theta = linspace(0,pi,numDir+1); theta(end) = [];
nu_list = [cos(theta); sin(theta)];
iIntN = findneighbours(r,numDir,nu_list,numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);

%% Finding neighbours for accurate scheme
[ind, axx, ayy, axy] = findneighboursaccurate(numBdy,numTot,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square,indicator_boundary);


%% Solving

% Matrices for the boundary conditions.
Mbc = speye(numBdy,numTot);
bc = uex(x(iBdy),y(iBdy));

% Matrices for Laplacian.
Dxx = SecondDerivMatrix(x,y,iInt,iIntN,[1;0],numDir,numInt,numTot,ones(numInt,1));
iIntNyy = findneighbours(r,1,[0;1],numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);
Dyy = SecondDerivMatrix(x,y,iInt,iIntNyy,[0;1],1,numInt,numTot,ones(numInt,1));

% Initial guess.
u0 = [Mbc;-(Dxx+Dyy)]\[bc;-sqrt(2*f(x(iInt),y(iInt)))];

% Solver parameters
tol = 1/N^2;
alphaMin = tol/N;
eps = 1/N^2;

% Monotone solver
uM = Newton(@(u)MongeAmpere(u,Mbc,bc,f(x(iInt),y(iInt)),numInt,iInt,iIntN,nu_list,numDir,numTot,numFr,x,y,eps),u0,tol,alphaMin);

% Filtered solver
% we start as initial guess with the solution of the monotone scheme
[Dxx,Dyy,Dxy] = SetupDerivsStandard(numBdy,numInt,numTot,ind,axx,axy,ayy);
epsilon = epsilon_choice(h,dTh); % epsilon - filtered scheme parameter
uF = Newton(@(u)MongeAmpereFiltered(u,Mbc,bc,f(x(iInt),y(iInt)),numInt,iInt,iIntN,nu_list,numDir,numTot,numFr,x,y,eps,Dxx,Dyy,Dxy,epsilon),uM,tol,alphaMin);
end
