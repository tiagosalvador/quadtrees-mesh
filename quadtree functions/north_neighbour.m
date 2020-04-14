function neighbour = north_neighbour(square,quadtree,vertices_square,x0,y0,x,y)
%
% This function determines the north neighbour of the square in quadtreed
%
% INPUTS:   square   Integer representing a square in the quadtree
%           quadtree Quadtree
%           (x0,y0)  Point on the north edge of the square
%
% OUTPUTS:  neighbour   North neighbour of square
% 
% there may be more than one north neighbour. Here, we assume there is only
% one

[root,pos] = parent_position(square,quadtree);

neighbour = 0;
if root
    switch pos
        case 1
            aux = north_neighbour(root,quadtree);
            if aux
                neighbour = (quadtree(aux,3));
                if not(neighbour)
                    neighbour = aux;
                end
            end
        case 2
            aux = north_neighbour(root,quadtree);
            if aux
                neighbour = (quadtree(aux,4));
                if not(neighbour)
                    neighbour = aux;
                end
            end
        case 3
            neighbour = (quadtree(root,1));
        otherwise
            neighbour = (quadtree(root,2));
    end
    
    
end

if and(neighbour,nargin > 2)
    while sum(quadtree(neighbour,[3 4])) > 0
        if x(vertices_square(quadtree(neighbour,4),4)) <= x0
            neighbour = quadtree(neighbour,4);
        else
            neighbour = quadtree(neighbour,3);
        end
    end  
end

end
