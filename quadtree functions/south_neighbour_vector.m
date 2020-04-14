function neighbour = south_neighbour_vector(square,quadtree,vertices_square,x0,y0,x,y)
% This function determines the south neighbour of the square in quadtreed
%
% INPUTS:   square   Integer representing a square in the quadtree
%           quadtree Quadtree
%
% OUTPUTS:  neighbour   South neighbour of square
% 
% there may be more than one south neighbour. Here, we assume there is only
% one

[root,pos] = parent_position(square,quadtree);

neighbour = zeros(size(square));

indices = and(root > 0,pos == 1);
neighbour(indices) = quadtree(root(indices),3);

indices = and(root > 0,pos == 2);
neighbour(indices) = quadtree(root(indices),4);

indices = and(root > 0,pos == 3);
if length(indices) > 0
    aux = south_neighbour_vector(root(indices),quadtree);
    indices = find(indices);
    if length(find(aux)) > 0
        aux2 = aux(aux>0);
        aux3 = quadtree(aux(aux>0),1);
        aux3(aux3==0) = aux2(aux3==0);
        neighbour(indices(aux>0)) = aux3;
    end
end

indices = and(root > 0,pos == 4);
if length(indices) > 0
    aux = south_neighbour_vector(root(indices),quadtree);
    indices = find(indices);
    if length(find(aux)) > 0
        aux2 = aux(aux>0);
        aux3 = quadtree(aux(aux>0),2);
        aux3(aux3==0) = aux2(aux3==0);
        neighbour(indices(aux>0)) = aux3;
    end
end

if nargin > 2
    
    indices = find(neighbour);
    aux =  sum(quadtree(neighbour(indices),[1 2]),2) > 0;
    while sum(aux) > 0
        aux2 = x(vertices_square(quadtree(neighbour(indices(aux)),1),2)) >= x0(indices(aux));
        indices_aux = indices(aux);
        neighbour(indices_aux(aux2)) = quadtree(neighbour(indices_aux(aux2)),1);
        neighbour(indices_aux(not(aux2))) = quadtree(neighbour(indices_aux(not(aux2))),2);
        indices = find(neighbour);
        aux =  sum(quadtree(neighbour(indices),[1 2]),2) > 0;
    end
end


end