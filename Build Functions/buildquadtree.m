function quadtree = buildquadtree(max_depth)

num_squares = (4^max_depth-1)/3; % total number of squares on quadtree

% The quadtree is structured in the following way: each row represents the
% following information of a square:
%       1st column - NW child square
%       2nd column - NE child square
%       3rd column - SW child square
%       4th column - SE child square
%       5th column - parent square

quadtree = zeros(num_squares,5);
quadtree(1,1:4) = [2 3 4 5];

quadtree(1:((4^(max_depth-1)-1)/3),1) = 2:4:num_squares;
quadtree(1:((4^(max_depth-1)-1)/3),2) = 3:4:num_squares;
quadtree(1:((4^(max_depth-1)-1)/3),3) = 4:4:num_squares;
quadtree(1:((4^(max_depth-1)-1)/3),4) = 5:4:num_squares;

aux = repmat(1:((4^(max_depth-1)-1)/3),[4 1]);
quadtree(2:((4^max_depth-1)/3),5) = aux(:);

end