function r = childrenQ(square,quadtree)
%
% This function determines if a square in quadtree has children
%
% INPUTS:   square   Integer representing a square in the quadtree
%           quadtree Quadtree
%
% OUTPUTS:  r

r = quadtree(square,1) > 0;

end
