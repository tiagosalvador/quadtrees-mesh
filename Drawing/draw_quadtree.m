function draw_quadtree(grid_parameters,mesh_parameters,insidedomain)

quadtree = mesh_parameters.quadtree;
vertices_square = mesh_parameters.vertices_square;
num_squares = size(quadtree,1);
xgrid = grid_parameters.xgrid;
ygrid = grid_parameters.ygrid;

for square = 1:num_squares
    if quadtree(square,1:4) == 0
    l = side_length(square,vertices_square,xgrid);
    vertice = vertices_square(square,4);
    aux = vertices_square(square,1:4);
    switch sum(insidedomain(xgrid(aux),ygrid(aux)))
        case 4
            rectangle_color = 'white';
            edge_color = 'black';
        case 0
            rectangle_color = 'white';
            edge_color = 'white';
        otherwise
            rectangle_color = [0.5 0.5 0.5];
            edge_color = 'black';
    end
    rectangle('Position', [xgrid(vertice) ygrid(vertice) l l],'FaceColor',rectangle_color,'EdgeColor',edge_color);
    end
end
end