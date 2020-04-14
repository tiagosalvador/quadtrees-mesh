function [grid_parameters,mesh_parameters] = buildsquares_boundary(dTh,grid_parameters,mesh_parameters,domain_function)
% algorithm discussed in the paper

xgrid = grid_parameters.xgrid;
ygrid = grid_parameters.ygrid;

indicator_boundary = @(x,y) domain_function(x,y);
insidedomain = @(x,y) domain_function(x,y)<0;

quadtree = mesh_parameters.quadtree;
vertices_square = mesh_parameters.vertices_square;
squares_vertex = mesh_parameters.squares_vertex;

% when building squares_boundary we are assuming that the boundary
% intersects each edge of a square with no children at most once


% The squares_boundary is structured in the following way: each row
% represents the following information of a square:
%       1st column - index in xB/x of the grid point at which the boundary
%       intersects the N edge
%       2nd column - index in xB/x of the grid point at which the boundary
%       intersects the E edge
%       3rd column - index in xB/x of the grid point at which the boundary
%       intersects the S edge
%       4th column - index in xB/x of the grid point at which the boundary
%       intersects the W edge

num_squares = size(quadtree,1);
squares_boundary = sparse(num_squares,5);
grid_vertices_boundary = find(indicator_boundary(xgrid,ygrid)==0);
numBdy = 0;
index_grid_vertices_boundary = (numBdy+1):(numBdy+length(grid_vertices_boundary));

numBdy = numBdy+length(grid_vertices_boundary);

for vertice_label = [1 2 3 4]
    switch vertice_label
        case 1
            edges_boundary = [2 3];
        case 2
            edges_boundary = [3 4];
        case 3
            edges_boundary = [1 2];
        case 4
            edges_boundary = [1 4];
    end
    current_square = squares_vertex(grid_vertices_boundary,vertice_label);
    
    for edge=edges_boundary
        index_aux = index_grid_vertices_boundary(current_square>0);
        xB(index_aux) = xgrid(grid_vertices_boundary(current_square>0));
        yB(index_aux) = ygrid(grid_vertices_boundary(current_square>0));
        xgrid2x(grid_vertices_boundary(current_square>0)) = index_grid_vertices_boundary(current_square>0);
        squares_boundary(current_square(current_square>0),edge) = index_grid_vertices_boundary(current_square>0);
    end   
end

squares = find(quadtree(:,1)==0);

aux = vertices_square(squares,1:4);
aux = sum(insidedomain(xgrid(aux),ygrid(aux)),2);
    
squares = squares(and(aux>0,aux<4));


for cardinal_direction = ['N' 'E' 'S' 'W']
    switch cardinal_direction
        case 'N'
            edge_vertices = [1 2];
        case 'E'
            edge_vertices = [2 3];
        case 'S'
            edge_vertices = [3 4];
        case 'W'
            edge_vertices = [4 1];
    end
    
    aux = vertices_square(squares,edge_vertices);
    %squares_to_update = squares(sum(insidedomain(xgrid(aux),ygrid(aux)),2)==1);
    %aux = aux(sum(insidedomain(xgrid(aux),ygrid(aux)),2)==1,:);
    squares_to_update = squares(prod(indicator_boundary(xgrid(aux),ygrid(aux)),2)<0);
    aux = aux(prod(indicator_boundary(xgrid(aux),ygrid(aux)),2)<0,:);
    
    switch cardinal_direction
        case 'N'
            adjacent_square = north_neighbour_vector(squares_to_update,quadtree);
            check = or(not(adjacent_square),not(squares_boundary(adjacent_square,3)));
        case 'E'
            adjacent_square = east_neighbour_vector(squares_to_update,quadtree);
            check = or(not(adjacent_square),not(squares_boundary(adjacent_square,4)));
        case 'S'
            adjacent_square = south_neighbour_vector(squares_to_update,quadtree);
            check = or(not(adjacent_square),not(squares_boundary(adjacent_square,1)));
        case 'W'
            adjacent_square = west_neighbour_vector(squares_to_update,quadtree);
            check = or(not(adjacent_square),not(squares_boundary(adjacent_square,2)));
    end
    if size(aux(check,:),1) > 1
        x_coordinate = mean(xgrid(aux(check,:)),2);
        y_coordinate = mean(ygrid(aux(check,:)),2);
    elseif size(aux(check,:),1) == 1
        x_coordinate = mean(xgrid(aux(check,:)));
        y_coordinate = mean(ygrid(aux(check,:)));
    else
        x_coordinate = xgrid(aux(check,1));
        y_coordinate = ygrid(aux(check,1));
    end
    switch cardinal_direction
        case {'N','S'}
            fun = @(x_min,x_max,y0) fzero(@(x) indicator_boundary(x,y0),[x_min x_max]);
            x_coordinate = arrayfun(fun,xgrid(aux(check,1)),xgrid(aux(check,2)),y_coordinate);
        case {'E','W'}
            fun = @(y_min,y_max,x0) fzero(@(y) indicator_boundary(x0,y),[y_min y_max]);
            y_coordinate = arrayfun(fun,ygrid(aux(check,1)),ygrid(aux(check,2)),x_coordinate);
    end
    
    % Guarantees no repetead boundary points
    [~,addBdy_select] = setdiff([x_coordinate y_coordinate],[xB' yB'],'rows','stable');
    if not(length(addBdy_select) == length(x_coordinate))
        [~,w,e] = intersect([x_coordinate y_coordinate],[xB' yB'],'rows','stable');
        x_coordinate = x_coordinate(addBdy_select);
        y_coordinate = y_coordinate(addBdy_select);
        numBdyaux = (numBdy+1):(numBdy+length(x_coordinate));
        xB(numBdyaux) = x_coordinate;
        yB(numBdyaux) = y_coordinate;
        numBdy = numBdy+length(numBdyaux);
        numBdyaux(addBdy_select) = numBdyaux;
        numBdyaux(w) = e;
    else
        numBdyaux = (numBdy+1):(numBdy+length(x_coordinate));
        xB(numBdyaux) = x_coordinate;
        yB(numBdyaux) = y_coordinate;
        numBdy = numBdy+length(numBdyaux);
        
    end
    
    
    switch cardinal_direction
        case 'N'
            squares_boundary(squares_to_update(check),1) = numBdyaux;
        case 'E'
            squares_boundary(squares_to_update(check),2) = numBdyaux;
        case 'S'
            squares_boundary(squares_to_update(check),3) = numBdyaux;
        case 'W'
            squares_boundary(squares_to_update(check),4) = numBdyaux;
    end
    switch cardinal_direction
        case 'N'
            indices_squares = find(not(check)).*and(adjacent_square(not(check)),squares_boundary(adjacent_square(not(check)),3));
            squares_boundary(squares_to_update(indices_squares),1) = squares_boundary(adjacent_square(indices_squares),3);
        case 'E'
            indices_squares = find(not(check)).*and(adjacent_square(not(check)),squares_boundary(adjacent_square(not(check)),4));
            squares_boundary(squares_to_update(indices_squares),2) = squares_boundary(adjacent_square(indices_squares),4);
        case 'S'
            indices_squares = find(not(check)).*and(adjacent_square(not(check)),squares_boundary(adjacent_square(not(check)),1));
            squares_boundary(squares_to_update(indices_squares),3) = squares_boundary(adjacent_square(indices_squares),1);
        case 'W'
            indices_squares = find(not(check)).*and(adjacent_square(not(check)),squares_boundary(adjacent_square(not(check)),2));
            squares_boundary(squares_to_update(indices_squares),4) = squares_boundary(adjacent_square(indices_squares),2);
    end
end


situation = zeros(size(squares));
vertice_label = zeros(size(squares));
vertice_labels_aux(:,1) = zeros(size(squares));
vertice_labels_aux(:,2) = zeros(size(squares));
min_angle = zeros(size(squares));
max_angle = zeros(size(squares));


% N edge

aux1 = squares_boundary(squares,1);
aux2 = squares_boundary(squares,2);
auxN = and(not(aux1==aux2),and(aux1,aux2));
vertice_label(auxN) = 2;
vertice_labels_aux(auxN,1) = 1;
vertice_labels_aux(auxN,2) = 2;
min_angle(auxN) = 0;
max_angle(auxN) = pi/2;
situation(auxN) = 1;

% E edge

aux1 = squares_boundary(squares,2);
aux2 = squares_boundary(squares,3);
auxE = and(not(auxN),and(not(aux1==aux2),and(aux1,aux2)));
vertice_label(auxE) = 3;
vertice_labels_aux(auxE,1) = 2;
vertice_labels_aux(auxE,2) = 3;
min_angle(auxE) = -pi/2;
max_angle(auxE) = 0;
situation(auxE) = 1;

% S edge
    
aux1 = squares_boundary(squares,3);
aux2 = squares_boundary(squares,4);
auxS = and(and(not(auxN),not(auxE)),and(not(aux1==aux2),and(aux1,aux2)));

vertice_label(auxS) = 4;
vertice_labels_aux(auxS,1) = 3;
vertice_labels_aux(auxS,2) = 4;
min_angle(auxS) = 0;
max_angle(auxS) = pi/2;
situation(auxS) = 1;

% W edge
aux1 = squares_boundary(squares,4);
aux2 = squares_boundary(squares,1);
auxW = and(and(and(not(auxN),not(auxE)),not(auxS)),and(not(aux1==aux2),and(aux1,aux2)));

vertice_label(auxW) = 1;
vertice_labels_aux(auxW,1) = 4;
vertice_labels_aux(auxW,2) = 1;
min_angle(auxW) = -pi/2;
max_angle(auxW) = 0;
situation(auxW) = 1;


aux = and(and(and(not(auxN),not(auxE)),not(auxS)),not(auxW));
aux = and(aux,arrayfun(@(square)length(setdiff(unique(squares_boundary(square,1:4)),0))>1,squares));
situation(aux) = 3;

indices = find(situation==1);

node(indices) = diag(vertices_square(squares(indices),vertice_label(indices)));
x0(indices) = xgrid(node(indices));
y0(indices) = ygrid(node(indices));


aux = indicator_boundary(x0(indices),y0(indices))>0;
situation(indices(aux)) = 2;

switch_vertice_label = repmat([3 4 1 2],length(vertice_label(indices(aux))),1);
vertice_label(indices(aux)) = diag(switch_vertice_label(:,vertice_label(indices(aux))))';

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
old = vertice_labels_aux(indices(aux),1);

vertice_labels_aux(indices(aux),1) = vertice_labels_aux(indices(aux),2);
vertice_labels_aux(indices(aux),2) = old;


node(indices(aux)) = diag(vertices_square(squares(indices(aux)),vertice_label(indices(aux))));
x0(indices(aux)) = xgrid(node(indices(aux)));
y0(indices(aux)) = ygrid(node(indices(aux)));

x_coordinates(:,1) = xB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),1))));
x_coordinates(:,2) = xB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),2))));

y_coordinates(:,1) = yB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),1))));
y_coordinates(:,2) = yB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),2))));



m_list = (y_coordinates-[y0(indices(aux));y0(indices(aux))]')./(x_coordinates-[x0(indices(aux));x0(indices(aux))]');
angles = atan(m_list);
min_angle(indices(aux)) = min(angles,[],2);
max_angle(indices(aux)) = max(angles,[],2);
if sum(not(insidedomain(x0(indices(aux)),y0(indices(aux)))))
    disp('oh boy')
end

clear x_coordinates y_coordinates


indices_to_update= indices(aux);
active = sum(squares_boundary(squares(indices(aux)),1:4)~=0,2)>2;
situation(indices_to_update(active)) = 3;



indices = find(situation==3);

aux = insidedomain(xgrid(vertices_square(squares(indices),:)),ygrid(vertices_square(squares(indices),:)));
index1 = aux(:,1);
index(index1) = 1;
index2 = and(not(index1),aux(:,2));
index(index2) = 2;
index3 = and(not(index1),and(not(index2),aux(:,3)));
index(index3) = 3;
index4 = and(not(index1),and(not(index2),and(not(index3),aux(:,4))));
index(index4) = 4;

x0(indices) = diag(xgrid(vertices_square(squares(indices),index)));
y0(indices) = diag(ygrid(vertices_square(squares(indices),index)));

aux = and(squares_boundary(squares(indices),1),squares_boundary(squares(indices),3));
vertice_labels_aux(indices(aux),1) = 3;
vertice_labels_aux(indices(aux),2) = 1;

x_coordinates = zeros(length(indices(aux)),2);
y_coordinates = zeros(length(indices(aux)),2);

if ~isempty(x_coordinates)
x_coordinates(:,1) = xB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),1))));
x_coordinates(:,2) = xB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),2))));

y_coordinates(:,1) = yB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),1))));
y_coordinates(:,2) = yB(diag(squares_boundary(squares(indices(aux)),vertice_labels_aux(indices(aux),2))));
end;


m_list = (y_coordinates-[y0(indices(aux));y0(indices(aux))]')./(x_coordinates-[x0(indices(aux));x0(indices(aux))]');
angles = atan(m_list);
min_angle(indices(aux)) = min(angles,[],2);
max_angle(indices(aux)) = max(angles,[],2);

aux2 = and(squares_boundary(squares(indices),2),squares_boundary(squares(indices),4));

aux = and(not(aux),aux2);
vertice_labels_aux(indices(aux),1) = 2;
vertice_labels_aux(indices(aux),2) = 4;
situation(indices(aux)) = 4;

min_angle(indices(and(aux,index'==1))) = -pi/2;
max_angle(indices((and(aux,index'==1)))) = atan((yB(squares_boundary(squares(indices((and(aux,index'==1)))),2))-y0(indices((and(aux,index'==1)))))./(xB(squares_boundary(squares(indices((and(aux,index'==1)))),2))-x0(indices((and(aux,index'==1))))));

min_angle(indices(and(aux,index'==3))) = -pi/2;
max_angle(indices((and(aux,index'==3)))) = atan((yB(squares_boundary(squares(indices((and(aux,index'==3)))),4))-y0(indices((and(aux,index'==3)))))./(xB(squares_boundary(squares(indices((and(aux,index'==3)))),4))-x0(indices((and(aux,index'==3))))));

clear angles

boundary_resolution = 50;

angles = arrayfun(@(a,b) linspace(a,b,boundary_resolution),min_angle,max_angle,'UniformOutput',false);
angles = cell2mat(angles);

active = situation > 0;

m_list = tan(angles);
x_coordinate = zeros(size(m_list));
y_coordinate = zeros(size(m_list));
for j=2:(size(m_list,2)-1)
    m = m_list(:,j);
    
    m = m(active);
    
    
    check = zeros(length(active),1);
    check(active) = abs(m)<1;
    x_min = zeros(size(squares(active)));
    x_max = zeros(size(squares(active)));
    y_min = zeros(size(squares(active)));
    y_max = zeros(size(squares(active)));
    
    x_min(and(active,check)) = xgrid(vertices_square(squares(and(active,check)),1));
    x_max(and(active,check)) = xgrid(vertices_square(squares(and(active,check)),2));
    
    check = zeros(length(active),1);
    check(active) = abs(m)>=1;
    
    y_min(and(active,check)) = ygrid(vertices_square(squares(and(active,check)),3));
    y_max(and(active,check)) = ygrid(vertices_square(squares(and(active,check)),2));
    x_min(and(active,check)) = (y_min(and(active,check))-y0(and(active,check))')./m(abs(m)>=1)+x0(and(active,check))';
    x_max(and(active,check)) = (y_max(and(active,check))-y0(and(active,check))')./m(abs(m)>=1)+x0(and(active,check))';
        
    fun = @(x0,y0,m,x_min,x_max) fzero(@(x) indicator_boundary(x,m*(x-x0)+y0),[x_min x_max]);
    x_coordinate(active,j) = arrayfun(fun,x0(active)',y0(active)',m,x_min(active),x_max(active));
    y_coordinate(active,j) = m.*(x_coordinate(active,j)-x0(active)')+y0(active)';
end


y_0 = arrayfun(@(a,b) linspace(a,b,boundary_resolution),ygrid(vertices_square(squares,3)),ygrid(vertices_square(squares,2)),'UniformOutput',false);
y_0 = cell2mat(y_0);
active = situation == 3;

for j=2:(size(y_0,2)-1)
    
    y_0aux = y_0(:,j);
    
    y_0aux = y_0aux(active);
    
    x_min = zeros(size(squares(active)));
    x_max = zeros(size(squares(active)));
    
    x_min(active) = xgrid(vertices_square(squares(active),1));
    x_max(active) = xgrid(vertices_square(squares(active),2));
        
    fun = @(y0,x_min,x_max) fzero(@(x) indicator_boundary(x,y0),[x_min x_max]);
    x_coordinate(active,j) = arrayfun(fun,y_0aux,x_min(active),x_max(active));
    y_coordinate(active,j) = y_0aux;
end


x_0 = arrayfun(@(a,b) linspace(a,b,boundary_resolution),xgrid(vertices_square(squares,1)),xgrid(vertices_square(squares,2)),'UniformOutput',false);
x_0 = cell2mat(x_0);
active = situation == 4;


for j=2:(size(x_0,2)-1)
    
    x_0aux = x_0(:,j);
    
    x_0aux = x_0aux(active);
    
    y_min = zeros(size(squares(active)));
    y_max = zeros(size(squares(active)));
    
    y_min(active) = ygrid(vertices_square(squares(active),3));
    y_max(active) = ygrid(vertices_square(squares(active),2));
        
    fun = @(x0,y_min,y_max) fzero(@(y) indicator_boundary(x0,y),[y_min y_max]);
    x_coordinate(active,j) = x_0aux;
    y_coordinate(active,j) = arrayfun(fun,x_0aux,y_min(active),y_max(active));
end

active = situation>0;

aux1 = (x_coordinate(active,2)'-xB(diag(squares_boundary(squares(active),vertice_labels_aux(active,1)')))).^2+(y_coordinate(active,2)'-yB(diag(squares_boundary(squares(active),vertice_labels_aux(active,1)')))).^2;
aux2 = (x_coordinate(active,2)'-xB(diag(squares_boundary(squares(active),vertice_labels_aux(active,2)')))).^2+(y_coordinate(active,2)'-yB(diag(squares_boundary(squares(active),vertice_labels_aux(active,2)')))).^2;

[~,choice] = min([aux1;aux2]);

x_coordinate(active,1) = xB(diag(squares_boundary(squares(active),diag(vertice_labels_aux(active,choice)))));
y_coordinate(active,1) = yB(diag(squares_boundary(squares(active),diag(vertice_labels_aux(active,choice)))));

x_coordinate(active,end) = xB(diag(squares_boundary(squares(active),diag(vertice_labels_aux(active,-choice+3)))));
y_coordinate(active,end) = yB(diag(squares_boundary(squares(active),diag(vertice_labels_aux(active,-choice+3)))));


active = situation > 0;

i = 1:(boundary_resolution-1);
arclength = cumsum(sqrt((x_coordinate(:,i+1)-x_coordinate(:,i)).^2+(y_coordinate(:,i+1)-y_coordinate(:,i)).^2),2);

arclength = [zeros(length(squares),1) arclength];


for i = 1:length(squares)
    if active(i)
        interior_mesh_points = find(insidedomain(xgrid',ygrid'));
        sz1 = length(interior_mesh_points);
        sz2 = length(x_coordinate(i,:));
        xInt = repmat(xgrid(interior_mesh_points),1,sz2);
        yInt = repmat(ygrid(interior_mesh_points),1,sz2);
        xBdy = repmat(x_coordinate(i,:),sz1,1);
        yBdy = repmat(y_coordinate(i,:),sz1,1);        
        delta = min(min(sqrt((xInt-xBdy).^2+(yInt-yBdy).^2)));
        
        r = linspace(arclength(i,1),arclength(i,end),ceil(arclength(i,end)/(2*2*delta*tan(dTh/2)))+2);
        x_aux = x_coordinate(i,(arrayfun(@(a) find(arclength(i,:)>=a,1),r)));
        y_aux = y_coordinate(i,(arrayfun(@(a) find(arclength(i,:)>=a,1),r)));
        x_aux(1) = [];
        x_aux(end) = [];
        y_aux(1) = [];
        y_aux(end) = [];
        
        % Guarantees no repetead boundary points
        addBdy = setdiff([x_aux' y_aux'],[xB' yB'],'rows');
        x_aux = addBdy(:,1);
        y_aux = addBdy(:,2);
        
        numBdyaux = (numBdy+1):(numBdy+length(x_aux));
        xB(numBdyaux) = x_aux;
        yB(numBdyaux) = y_aux;
        numBdy = numBdy+length(numBdyaux);
        
        squares_boundary(squares(i), (4+1):(4+length(x_aux))) = numBdyaux;
    end
end

depth_squares = depth_level(squares,quadtree);
depth = max(depth_squares);
go = 1;
while depth >= 1
    clear indices
    squares_to_update = squares(depth_squares==depth);
    [parents, pos] = parent_position(squares_to_update,quadtree);
    squares_to_update = squares_to_update(parents>0);
    pos = pos(parents>0);
    parents = parents(parents>0);
    
    squares_boundary(parents(pos==1),1) = max(squares_boundary(parents(pos==1),1),squares_boundary(squares_to_update(pos==1),1));
    squares_boundary(parents(pos==1),4) = max(squares_boundary(parents(pos==1),4),squares_boundary(squares_to_update(pos==1),4));
    indices(pos==1,1) = 2;
    indices(pos==1,2) = 3;
    
    squares_boundary(parents(pos==2),1) = max(squares_boundary(parents(pos==2),1),squares_boundary(squares_to_update(pos==2),1));
    squares_boundary(parents(pos==2),2) = max(squares_boundary(parents(pos==2),2),squares_boundary(squares_to_update(pos==2),2));
    indices(pos==2,1) = 3;
    indices(pos==2,2) = 4;
    
    squares_boundary(parents(pos==3),3) = max(squares_boundary(parents(pos==3),3),squares_boundary(squares_to_update(pos==3),3));
    squares_boundary(parents(pos==3),4) = max(squares_boundary(parents(pos==3),4),squares_boundary(squares_to_update(pos==3),4));
    indices(pos==3,1) = 1;
    indices(pos==3,2) = 2;
    
    squares_boundary(parents(pos==4),2) = max(squares_boundary(parents(pos==4),2),squares_boundary(squares_to_update(pos==4),2));
    squares_boundary(parents(pos==4),3) = max(squares_boundary(parents(pos==4),3),squares_boundary(squares_to_update(pos==4),3));
    indices(pos==4,1) = 1;
    indices(pos==4,2) = 4;
    
    unique_parents = unique(parents);
    for k=1:length(unique_parents)
        
        parent = unique_parents(k);
        
        indices_parent = find(parents == parent);
        aux1 = sub2ind(size(squares_boundary),squares_to_update(indices_parent),indices(indices_parent,1));
        aux2 = sub2ind(size(squares_boundary),squares_to_update(indices_parent),indices(indices_parent,2));
        aux1 = union(aux1,aux2);
        
        I = repmat(squares_to_update(indices_parent),1,size(squares_boundary,2)-5+1);
        J = repmat(5:size(squares_boundary,2),length(squares_to_update(indices_parent)),1);
        
        aux2 = sub2ind(size(squares_boundary),I,J);
        aux = union(aux1(:),aux2(:));
        aux = unique(nonzeros(squares_boundary(aux)));
        squares_boundary(parent,5:(5+length(aux)-1)) = aux;
    end
    squares = union(unique_parents,setdiff(squares,squares_to_update));
    depth_squares = depth_level(squares,quadtree);
    depth = max(depth_squares);
end
mesh_parameters.squares_boundary = squares_boundary;

grid_parameters.numBdy = numBdy;
grid_parameters.xgrid2x = xgrid2x;
grid_parameters.xB = xB;
grid_parameters.yB = yB;
end

