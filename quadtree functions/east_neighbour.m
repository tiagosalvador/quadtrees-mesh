function neighbour = east_neighbour(square,quadtree,vertices_square,x0,y0,x,y)
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

neighbour = 0;
if root
    switch pos
        case 1
            neighbour = quadtree(root,2);
        case 2
            aux = east_neighbour(root,quadtree);
            if aux
                neighbour = quadtree(aux,1);
                if not(neighbour)
                    neighbour = aux;
                end
            end
        case 3
            neighbour = quadtree(root,4);
        otherwise
            aux = east_neighbour(root,quadtree);
            if aux
                neighbour = quadtree(aux,3);
                if not(neighbour)
                    neighbour = aux;
                end
            end
    end
end
if and(neighbour,nargin > 2)
    while sum(quadtree(neighbour,[1 3])) > 0
        if y(vertices_square(quadtree(neighbour,3),1)) >= y0
            neighbour = quadtree(neighbour,3);
        else
            neighbour = quadtree(neighbour,1);
        end
    end  
end
end
