function [uM,uA,grid_parameters,mesh_parameters,gridR_parameters,meshR_parameters] = solveConvEnvAdaptivity(a,b,domain_function,max_depth,uex,g,dTh_choice)
% Solves the Convex Envelopeequation on [a,b]x[a,b] on a uniform grid of
% size n by n.
% Args:
%    domain_function: variable to determine the domain
%    max_depth: maximum depth of the initial quadtree
%    uex: exact solution to impose boundary condtions
%    g: obstacle function in the Convex Envelope equation
%    dTh_choice: function to determine dTh
%    epsilon_choice: function to determine epsilon for the filtered scheme
% Returns:
%    uM: approximate solution on stardard grid
%    uA: approximate solution on adaptive grid
%    grid_paramters: Matlab structure containing grid information
%    mesh_paramters: Matlab structure containing mesh information
%    gridR_paramters: Matlab structure containing refined grid information
%    meshR_paramters: Matlab structure containing refined mesh information

%% Setup

insidedomain = @(x,y) domain_function(x,y)<0;

[grid_parameters,dTh,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice);
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

%% finding neighbours for monotone scheme
r = h*(1+cos(dTh/2)*cot(dTh/2)+sin(dTh/2));
numDir = round(pi/dTh);
theta = linspace(0,pi,numDir+1); theta(end) = [];
nu_list = [cos(theta); sin(theta)];
iIntN = findneighbours(r,numDir,nu_list,numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);

%% Solving

% Matrices for the boundary conditions.
Mbc = speye(numBdy,numTot);
bc = uex(x(iBdy),y(iBdy));

% Initial guess
u = g(x,y);

% Solve.
tol = 1e-8;
u = CEActive(u,Mbc,bc,g(x(iInt),y(iInt)),numInt,iInt,iBdy,iIntN,nu_list,numDir,numTot,x,y,tol);
uM = u;


%% Refinement

% Determine Dxx,Dxy,Dyy
Dxx = SecondDerivMatrix(x,y,iInt,iIntN,[1;0],numDir,numInt,numTot,ones(numInt,1));
iIntNyy = findneighbours(r,1,[0;1],numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);
Dyy = SecondDerivMatrix(x,y,iInt,iIntNyy,[0;1],1,numInt,numTot,ones(numInt,1));
iIntNxy = findneighbours(r,1,[cos(pi/4);sin(pi/4)],numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);
Dxy = SecondDerivMatrix(x,y,iInt,iIntNxy,[0;1],1,numInt,numTot,ones(numInt,1));

% Refine quadtree
[gridR_parameters,meshR_parameters] = ...
    refine_balance_quadtree_posteriori(grid_parameters,mesh_parameters,domain_function,Dxx,Dyy,Dxy,u);

%% Building squares_boundary
[gridR_parameters,meshR_parameters] = buildsquares_boundary(dTh,gridR_parameters,meshR_parameters,domain_function);

%% Grid parameters
xgrid = gridR_parameters.xgrid;
ygrid = gridR_parameters.ygrid;
xgrid2x = gridR_parameters.xgrid2x;
numTotGrid = gridR_parameters.numTotGrid;
numBdy = gridR_parameters.numBdy;
xB = gridR_parameters.xB;
yB = gridR_parameters.yB;

in = insidedomain(xgrid,ygrid);
xInt = xgrid(in);
yInt = ygrid(in);

numInt = length(xInt(:));
x = [xB(:);xInt(:)];
y = [yB(:);yInt(:)];
iBdy = (1:numBdy)';
aux = find(in);
xgrid2x(aux) = numBdy + (1:length(aux));
xgrid2x((end+1):(end+(numTotGrid-length(xgrid2x)))) = 0;
numTot = numBdy + numInt;
iInt = (numBdy+1:numTot)';


gridR_parameters.xgrid2x = xgrid2x;
gridR_parameters.iBdy = iBdy;
gridR_parameters.iInt = iInt;
gridR_parameters.x = x;
gridR_parameters.y = y;
gridR_parameters.numTotGrid = numTotGrid;
gridR_parameters.numInt = numInt;
gridR_parameters.numTot = numTot;

%% Finding neighbours for refined scheme

quadtree = meshR_parameters.quadtree;
squares_boundary = meshR_parameters.squares_boundary;
squares_vertex = meshR_parameters.squares_vertex;
vertices_square = meshR_parameters.vertices_square;

lengthNW = side_length(squares_vertex(in,1),vertices_square,xgrid);
lengthNE = side_length(squares_vertex(in,2),vertices_square,xgrid);
lengthSW = side_length(squares_vertex(in,3),vertices_square,xgrid);
lengthSE = side_length(squares_vertex(in,4),vertices_square,xgrid);

locallength = max(lengthNW,max(lengthNE,max(lengthSW,lengthSE)));

[hrefined,~,refinement_level] = unique(locallength);

%% Finding neighbours for monotone scheme in refined grid
dTh = dTh_choice(hrefined);
r = hrefined.*(1+cos(dTh/2).*cot(dTh/2)+sin(dTh/2));
numDir = round(pi./dTh);

clear nu_list

theta = arrayfun(@(numDir) linspace(0,pi,numDir+1),numDir,'UniformOutput',false);

for kk = 1:length(hrefined)
    theta{kk}(end) = [];
    nu_list{kk} = [cos(theta{kk}); sin(theta{kk})];
end
iIntN = findneighboursrefined(refinement_level,r,numDir,nu_list,numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);

%% Solving

% Matrices for the boundary conditions.
Mbc = speye(numBdy,numTot);
bc = uex(x(iBdy),y(iBdy));

% Initial guess
u = g(x,y);

% Solve.
tol = 1e-8;
u = CEActiverefined(refinement_level,u,Mbc,bc,g(x(iInt),y(iInt)),numInt,iInt,iBdy,iIntN,nu_list,numDir,numTot,x,y,tol);
uA = u;

end