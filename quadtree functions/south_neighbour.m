function neighbour = south_neighbour(square,quadtree,vertices_square,x0,y0,x,y)
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

neighbour = 0;
if root
    switch pos
        case 1
            neighbour = quadtree(root,3);
        case 2
            neighbour = quadtree(root,4);
        case 3
            aux = south_neighbour(root,quadtree);
            if aux
                neighbour = quadtree(aux,1);
                if not(neighbour)
                    neighbour = aux;
                end
            end
        otherwise
            aux = south_neighbour(root,quadtree);
            if aux
                neighbour = quadtree(aux,2);
                if not(neighbour)
                    neighbour = aux;
                end
            end
    end
end
if and(neighbour,nargin > 2)
    while sum(quadtree(neighbour,[1 2])) > 0
        if x(vertices_square(quadtree(neighbour,1),2)) >= x0
            neighbour = quadtree(neighbour,1);
        else
            neighbour = quadtree(neighbour,2);
        end
    end  
end
end
