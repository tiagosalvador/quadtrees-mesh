function [quadtree,num_squares,vertices_square,...
    squares_vertex,xgrid,ygrid,numTotGrid] = ...
    refine_balance_quadtree(quadtree,num_squares,vertices_square,squares_vertex,xgrid,ygrid,numTotGrid,indicator_boundary,insidedomain,refinement_function)

h = ygrid(2)-ygrid(1);
squares_to_refine = find(quadtree(:,1)==0);
aux = vertices_square(squares_to_refine,1:4);
%squares_to_refine = squares_to_refine(and(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>0,sum(insidedomain(xgrid(aux),ygrid(aux)),2)<4));
squares_to_refine = squares_to_refine(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>0);
aux = vertices_square(squares_to_refine,1:4);
squares_to_refine = squares_to_refine(sum(abs(refinement_function(xgrid(aux),ygrid(aux)))<0.1,2)>0);
k = 0;
while k<2
    for current_square = squares_to_refine'
        [quadtree,num_squares,vertices_square,...
            squares_vertex,xgrid,ygrid,numTotGrid] = ...
            refine_square(current_square,quadtree,num_squares,vertices_square,squares_vertex,xgrid,ygrid,numTotGrid);
    end
    squares_to_refine = quadtree(squares_to_refine,1:4);
    squares_to_refine = squares_to_refine(:);
    k=k+1;
end

%% balancing quadtree

squares_to_refine = find(quadtree(:,1)==0);
aux = vertices_square(squares_to_refine,1:4);
squares_to_refine = squares_to_refine(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>0);

while squares_to_refine
    current_square = squares_to_refine(1);
    depth_current_square = depth_level(current_square,quadtree);
    %squares_to_refine(1) = [];
    squares_to_refine = setdiff(squares_to_refine,current_square);
    refineQ = 0;
    for cardinal_direction = ['N' 'E' 'S' 'W']
        switch cardinal_direction
            case 'N'
                neighbour = north_neighbour(current_square,quadtree);
                location1 = 3;
                location2 = 4;
            case 'E'
                neighbour = east_neighbour(current_square,quadtree);
                location1 = 1;
                location2 = 3;
            case 'S'
                neighbour = south_neighbour(current_square,quadtree);
                location1 = 1;
                location2 = 2;
            case 'W'
                neighbour = west_neighbour(current_square,quadtree);
                location1 = 2;
                location2 = 4;
                
        end
        if neighbour
            aux = vertices_square(neighbour,1:4);
            if sum(insidedomain(xgrid(aux),ygrid(aux)))>0
                depth_neighbour = depth_level(neighbour,quadtree);
                if depth_neighbour >= depth_current_square
                    if childrenQ(neighbour,quadtree)
                        if or(childrenQ(quadtree(neighbour,location1),quadtree),childrenQ(quadtree(neighbour,location2),quadtree))
                            refineQ = 1;
                        end
                    end
                end
            end
        end
        if refineQ
            break
        end
    end
    if refineQ
        [quadtree,num_squares,vertices_square,...
            squares_vertex,xgrid,ygrid,numTotGrid] = ...
            refine_square(current_square,quadtree,num_squares,vertices_square,squares_vertex,xgrid,ygrid,numTotGrid);
        squares_to_refine = [squares_to_refine;quadtree(current_square,1:4)'];
        for cardinal_direction = ['N' 'E' 'S' 'W']
            switch cardinal_direction
                case 'N'
                    neighbour = north_neighbour(current_square,quadtree);
                case 'E'
                    neighbour = east_neighbour(current_square,quadtree);
                case 'S'
                    neighbour = south_neighbour(current_square,quadtree);

                case 'W'
                    neighbour = west_neighbour(current_square,quadtree);                    
            end
            if neighbour
                if not(childrenQ(neighbour,quadtree))
                    aux = vertices_square(neighbour,1:4);
                    if sum(insidedomain(xgrid(aux),ygrid(aux)))>0
                        squares_to_refine = [squares_to_refine;neighbour];
                    end
                end
            end
        end
        
    end
    
end

% This would ensure that every neighbour square of a boundary square was at
% most at the same depth
% go = 1;
% while go
%     go = 0;
%     squares_to_refine = [];
%     
%     for cardinal_direction = ['N' 'E' 'S' 'W']
%         squares = find(quadtree(:,1)==0);
%         aux = vertices_square(squares,1:4);
%         squares = squares(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>3);
%         
%         switch cardinal_direction
%             case 'N'
%                 neighbours = north_neighbour_vector(squares,quadtree);
%                 
%             case 'E'
%                 neighbours = east_neighbour_vector(squares,quadtree);
%                 
%             case 'S'
%                 neighbours = south_neighbour_vector(squares,quadtree);
%                 
%             case 'W'
%                 neighbours = west_neighbour_vector(squares,quadtree);
%         end
%         aux = vertices_square(neighbours,1:4);
%         neighbours = neighbours(and(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>0,sum(insidedomain(xgrid(aux),ygrid(aux)),2)<4));
%         squares = squares(and(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>0,sum(insidedomain(xgrid(aux),ygrid(aux)),2)<4));
%         
%         depth_squares = depth_level(squares,quadtree);
%         depth_neighbours = depth_level(neighbours,quadtree);
%         
%         squares_to_refine = [squares_to_refine;neighbours(depth_neighbours<depth_squares)];
%         squares_to_refine = unique(squares_to_refine);
%     end
%     
%     for current_square = squares_to_refine'
%         [quadtree,num_squares,vertices_square,...
%             squares_vertex,xgrid,ygrid,numTotGrid] = ...
%             refine_square(current_square,quadtree,num_squares,vertices_square,squares_vertex,xgrid,ygrid,numTotGrid);
%     end
%     if not(isempty(squares_to_refine))
%         go = 1;
%     end
%     % we might want to balance the quadtree as well
% end
end