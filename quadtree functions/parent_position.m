function [root,pos] = parent_position(square,quadtree)
%
% This function determines the parent of the square in quadtreed
%
% INPUTS:   square   Integer representing a square in the quadtree
%           quadtree Quadtree
%
% OUTPUTS:  root     Parent square
%           pos      Position of square in the parent:
%                    1 = NW, 2 = NE, 3 = SE, 4 = SW      
%

root = quadtree(square,5);
pos = zeros(size(square));

if isempty(root)
    root = [];
end;
if isempty(pos)
    pos = [];
end;

if ~isempty(root(root>0))
    [I,J] = find(quadtree(root(root>0),:)==repmat(square(root>0),1,5));
    
    aux(I) = J;
    pos(root>0) = aux';
end

%pos(root>0) = arrayfun(@(x,y)find(quadtree(x,:) == y),root(root>0),square(root > 0));

% if root
%     pos = find(quadtree(root,:) == square);        
% else
%     pos = 0;
% end
end

