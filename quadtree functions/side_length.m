function l = side_length(square,vertices_square,xgrid)
l = xgrid(vertices_square(square,2))-xgrid(vertices_square(square,1));
end