function paint_square(current_square,vertices_square,xgrid,ygrid)
l = side_length(current_square,vertices_square,xgrid);
vertice = vertices_square(current_square,4);
rectangle('Position', [xgrid(vertice) ygrid(vertice) l l],'FaceColor','blue');
end