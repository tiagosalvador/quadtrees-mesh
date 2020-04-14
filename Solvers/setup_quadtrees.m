function [grid_parameters,dTh,mesh_parameters] = setup_quadtrees(a,b,domain_function,max_depth,dTh_choice)

insidedomain = @(x,y) domain_function(x,y)<0;

grid_parameters = struct;

N = 2^(max_depth-1)+1; % total number of points on the grid
h = (b-a)/(N-1); % length of the smallest squares on the initial quadtree
[xgrid,ygrid] = meshgrid(linspace(a,b,N)); % grid
xgrid = xgrid(:);
ygrid = ygrid(:);
grid_parameters.xgrid = xgrid;
grid_parameters.ygrid = ygrid;
numTotGrid = N^2; % total number of grid points
grid_parameters.numTotGrid = numTotGrid;


dTh = dTh_choice(h); % dTheta - montone scheme parameter

%% Building quadtree
quadtree = buildquadtree(max_depth);

%% Building vertices_square & squares_vertex
[vertices_square, squares_vertex] = buildvertices_square_and_squares_vertex(a,b,h,max_depth,N,quadtree,xgrid,ygrid);

mesh_parameters.quadtree = quadtree;
mesh_parameters.vertices_square = vertices_square;
mesh_parameters.squares_vertex = squares_vertex;

%% Building squares_boundary
[grid_parameters,mesh_parameters] = buildsquares_boundary(dTh,grid_parameters,mesh_parameters,domain_function);

%% Grid parameters
xgrid = grid_parameters.xgrid;
ygrid = grid_parameters.ygrid;
xgrid2x = grid_parameters.xgrid2x;
numTotGrid = grid_parameters.numTotGrid;
numBdy = grid_parameters.numBdy;
xB = grid_parameters.xB;
yB = grid_parameters.yB;

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

grid_parameters.N = N;
grid_parameters.h = h;
grid_parameters.xgrid2x = xgrid2x;
grid_parameters.iBdy = iBdy;
grid_parameters.iInt = iInt;
grid_parameters.x = x;
grid_parameters.y = y;
grid_parameters.numInt = numInt;
grid_parameters.numBdy = numBdy;
grid_parameters.numTot = numTot;

end