function neighbour = west_neighbour(square,quadtree,vertices_square,x0,y0,x,y)
%
% This function determines the west neighbour of the square in quadtreed
%
% INPUTS:   square   Integer representing a square in the quadtree
%           quadtree Quadtree
%
% OUTPUTS:  neighbour   west neighbour of square
% 
% there may be more than one west neighbour. Here, we assume there is only
% one

[root,pos] = parent_position(square,quadtree);

neighbour = 0;
if root
    switch pos
        case 1
            aux = west_neighbour(root,quadtree);
            if aux
                neighbour = quadtree(aux,2);
                if not(neighbour)
                    neighbour = aux;
                end
            end
        case 2
            neighbour = quadtree(root,1);
        case 3
            aux = west_neighbour(root,quadtree);
            if aux
                neighbour = quadtree(aux,4);
                if not(neighbour)
                    neighbour = aux;
                end
            end
            
        case 4
            neighbour = quadtree(root,3);
    end
end
if and(neighbour,nargin > 2)
    while sum(quadtree(neighbour,[2 4])) > 0
        if y(vertices_square(quadtree(neighbour,2),3)) <= y0
            neighbour = quadtree(neighbour,2);
        else
            neighbour = quadtree(neighbour,4);
        end
    end
end
end
