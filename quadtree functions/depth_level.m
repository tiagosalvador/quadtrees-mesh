function depth = depth_level(squares,quadtree)
%
% This function determines the depth level at which a square is
%
% INPUTS:  square   Integer representing a square in the quadtree
%          quadtree Quadtree
% 
% OUTPUTS: depth    Depth level of the square in the quadtree
depth = zeros(size(squares));
current_squares = quadtree(squares,5);
while sum(current_squares)    
    indices = current_squares>0;
    current_squares(indices) = quadtree(current_squares(indices),5);
    depth(indices) = depth(indices)+1;
end


end