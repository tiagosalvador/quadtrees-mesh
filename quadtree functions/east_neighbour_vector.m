function neighbour = east_neighbour_vector(square,quadtree,vertices_square,x0,y0,x,y)
%
% This function determines the north neighbour of the square in quadtreed
%
% INPUTS:   square   Integer representing a square in the quadtree
%           quadtree Quadtree
%
% OUTPUTS:  neighbour   east neighbour of square
% 
% there may be more than one east neighbour. Here, we assume there is only
% one

[root,pos] = parent_position(square,quadtree);


neighbour = zeros(size(square));

indices = and(root > 0,pos == 1);
neighbour(indices) = quadtree(root(indices),2);

indices = and(root > 0,pos == 2);
if length(indices) > 0
    aux = east_neighbour_vector(root(indices),quadtree);
    indices = find(indices);
    if length(find(aux)) > 0
        aux2 = aux(aux>0);
        aux3 = quadtree(aux(aux>0),1);
        aux3(aux3==0) = aux2(aux3==0);
        neighbour(indices(aux>0)) = aux3;
    end
end
indices = and(root > 0,pos == 3);
neighbour(indices) = quadtree(root(indices),4);

indices = and(root > 0,pos == 4);
if length(indices) > 0
    aux = east_neighbour_vector(root(indices),quadtree);
    indices = find(indices);
    if length(find(aux)) > 0
        aux2 = aux(aux>0);
        aux3 = quadtree(aux(aux>0),3);
        aux3(aux3==0) = aux2(aux3==0);
        neighbour(indices(aux>0)) = aux3;
    end
end

if nargin > 2
    indices = find(neighbour);
    aux =  sum(quadtree(neighbour(indices),[1 3]),2) > 0;
    while sum(aux) > 0
        aux2 = y(vertices_square(quadtree(neighbour(indices(aux)),3),1)) >= y0(indices(aux));
        indices_aux = indices(aux);
        neighbour(indices_aux(aux2)) = quadtree(neighbour(indices_aux(aux2)),3);
        neighbour(indices_aux(not(aux2))) = quadtree(neighbour(indices_aux(not(aux2))),1);
        indices = find(neighbour);
        aux =  sum(quadtree(neighbour(indices),[1 3]),2) > 0;
    end
    
end


end
