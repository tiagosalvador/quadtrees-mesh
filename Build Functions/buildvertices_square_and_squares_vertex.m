function [vertices_square, squares_vertex] = ...
    buildvertices_square_and_squares_vertex(a,b,h,max_depth,N,quadtree,xgrid,ygrid)


% The vertices_square is structured in the following way: each row
% represents the following information of a square:
%       1st column - index in xgrid/ygrid of the NW vertex
%       2nd column - index in xgrid/ygrid of the NE vertex
%       3rd column - index in xgrid/ygrid of the SE vertex
%       4th column - index in xgrid/ygrid of the SW vertex

% The squares_vertex is structured in the following way: each row
% represents the following information of a vertex (index in xgrid/ygrid):
%       1st column - NW square
%       2nd column - NE square
%       3rd column - SW square
%       4th column - SE square

num_squares = size(quadtree,1);
vertices_square = zeros(num_squares,4);
squares_vertex = zeros(N^2,4);

vertices_square(1,1:4) = [N N^2 N^2-N+1 1];

squares_vertex(vertices_square(1,1),4) = 1;
squares_vertex(vertices_square(1,2),3) = 1;
squares_vertex(vertices_square(1,3),1) = 1;
squares_vertex(vertices_square(1,4),2) = 1;



squares = (2:num_squares)';
depth_squares = depth_level(squares,quadtree);

depth = 1;
while depth<= max_depth-1
    current_square = squares(depth_squares==depth);
    [root,pos] = parent_position(current_square,quadtree);
    dx = (xgrid(vertices_square(root,2))-xgrid(vertices_square(root,1)))/2;
    aux_vertice = [xgrid(vertices_square(root,1)) ygrid(vertices_square(root,1))];
    
    aux_vertice(pos==1,1) = aux_vertice(pos==1,1);
    aux_vertice(pos==1,2) = aux_vertice(pos==1,2) - dx(pos==1);
    
    aux_vertice(pos==2,1) = aux_vertice(pos==2,1)+dx(pos==2);
    aux_vertice(pos==2,2) = aux_vertice(pos==2,2)-dx(pos==2);
    
    aux_vertice(pos==3,1) = aux_vertice(pos==3,1);
    aux_vertice(pos==3,2) = aux_vertice(pos==3,2)-2*dx(pos==3);
    
    aux_vertice(pos==4,1) = aux_vertice(pos==4,1)+dx(pos==4);
    aux_vertice(pos==4,2) = aux_vertice(pos==4,2)-2*dx(pos==4);
    
    % vertice 1 - NW
    x_coordinate = aux_vertice(:,1);
    y_coordinate = aux_vertice(:,2)+dx;
    vertices_square(current_square,1) = sub2ind([N N],round((y_coordinate-a)/h+1),round((x_coordinate+0-a)/h+1));
    
    % vertice 2 - NE
    x_coordinate = aux_vertice(:,1)+dx;
    y_coordinate = aux_vertice(:,2)+dx;
    vertices_square(current_square,2) = sub2ind([N N],round((y_coordinate-a)/h+1),round((x_coordinate+0-a)/h+1));
    
    % vertice 3 - SE
    x_coordinate = aux_vertice(:,1)+dx;
    y_coordinate = aux_vertice(:,2);
    vertices_square(current_square,3) = sub2ind([N N],round((y_coordinate-a)/h+1),round((x_coordinate+0-a)/h+1));
    
    % vertice 4 - SW
    x_coordinate = aux_vertice(:,1);
    y_coordinate = aux_vertice(:,2);
    vertices_square(current_square,4) = sub2ind([N N],round((y_coordinate-a)/h+1),round((x_coordinate+0-a)/h+1));
    
    squares_vertex(vertices_square(current_square,1),4) = current_square;
    squares_vertex(vertices_square(current_square,2),3) = current_square;
    squares_vertex(vertices_square(current_square,3),1) = current_square;
    squares_vertex(vertices_square(current_square,4),2) = current_square;
    
    depth = depth+1;
end

end