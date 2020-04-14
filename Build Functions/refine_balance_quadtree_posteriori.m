function [gridR_parameters,meshR_parameters] = ...
    refine_balance_quadtree_posteriori(grid_parameters,mesh_parameters,domain_function,Dxx,Dyy,Dxy,u);

insidedomain = @(x,y) domain_function(x,y)<0;

gridR_parameters = grid_parameters;
meshR_parameters = mesh_parameters;

quadtree = mesh_parameters.quadtree;
num_squares = size(quadtree,1);
vertices_square = mesh_parameters.vertices_square;
squares_vertex = mesh_parameters.squares_vertex;

xgrid = grid_parameters.xgrid;
ygrid = grid_parameters.ygrid;
numTotGrid = grid_parameters.numTotGrid;
xgrid2x = grid_parameters.xgrid2x;


[~,J,~] = find(xgrid2x);
in = insidedomain(xgrid(J),ygrid(J));
uxxgrid = zeros(size(xgrid));
uxxgrid(J(in)) = Dxx*u;
uyygrid = zeros(size(xgrid));
uyygrid(J(in)) = Dyy*u;
uxygrid = zeros(size(xgrid));
uxygrid(J(in)) = Dxy*u;
normhessiangrid = sqrt(uxxgrid.^2+uyygrid.^2+2*uxygrid.^2);
h = ygrid(2)-ygrid(1);
squares_to_refine = find(quadtree(:,1)==0);
aux = vertices_square(squares_to_refine,1:4);
squares_to_refine = squares_to_refine(sum(insidedomain(xgrid(aux),ygrid(aux)),2)>0);
aux = vertices_square(squares_to_refine,1:4);
squares_to_refine = squares_to_refine(sum(h*normhessiangrid(aux)>0.5,2)>0);
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

gridR_parameters.xgrid = xgrid;
gridR_parameters.ygrid = ygrid;
gridR_parameters.numTotGrid = numTotGrid;

meshR_parameters.quadtree = quadtree;
meshR_parameters.vertices_square = vertices_square;
meshR_parameters.squares_vertex = squares_vertex;

end