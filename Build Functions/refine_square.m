function [quadtree,num_squares,vertices_square,...
    squares_vertex,xgrid,ygrid,numTotGrid] = ...
    refine_square(current_square,quadtree,num_squares,vertices_square,squares_vertex,xgrid,ygrid,numTotGrid)

square_NW = num_squares+1;
square_NE = num_squares+2;
square_SW = num_squares+3;
square_SE = num_squares+4;

squares_aux = [square_NW square_NE square_SW square_SE];

quadtree(current_square,1:4) = (num_squares+1):(num_squares+4);
quadtree((num_squares+1):(num_squares+4),5) = current_square;
num_squares = num_squares+4;



index_NW = vertices_square(current_square,1);
index_NE = vertices_square(current_square,2);
index_SE = vertices_square(current_square,3);
index_SW = vertices_square(current_square,4);



squares_vertex(index_NW,4) = square_NW;
squares_vertex(index_NE,3) = square_NE;
squares_vertex(index_SW,2) = square_SW;
squares_vertex(index_SE,1) = square_SE;

for cardinal_direction = ['N' 'E' 'S' 'W']
    switch cardinal_direction
        case 'N'
            neighbour = north_neighbour(current_square,quadtree);
            index1 = 1; % where direction1 is stored
            index2 = 2; % where direction2 is stored
            direction1 = 3; % SW
            direction2 = 4; % SE
            index = 3;
        case 'E'
            neighbour = east_neighbour(current_square,quadtree);
            index1 = 2;
            index2 = 4;
            direction1 = 1; % NW
            direction2 = 3; % SW
            index = 4;
        case 'S'
            neighbour = south_neighbour(current_square,quadtree);
            index1 = 3;
            index2 = 4;
            direction1 = 1; % NW
            direction2 = 2; % NE
            index = 2;
        case 'W'
            neighbour = west_neighbour(current_square,quadtree);
            index1 = 1;
            index2 = 3;
            direction1 = 2; % NE
            direction2 = 4; % SE
            index = 3;
    end
    
    if neighbour
        if childrenQ(neighbour,quadtree)
            children = quadtree(neighbour,direction1);
            index = vertices_square(children,index);
            switch cardinal_direction
                case {'N','S'}
                    xmin = min(xgrid(vertices_square(children,direction1)),xgrid(vertices_square(children,direction2)));
                    xmax = max(xgrid(vertices_square(children,direction1)),xgrid(vertices_square(children,direction2)));
                    indices = find(and(ygrid==ygrid(index),and(xgrid>xmin,xgrid<xmax)));
                    squares_vertex(indices,direction1) = squares_aux(index1);
                    squares_vertex(indices,direction2) = squares_aux(index1);
                case {'W','E'}
                    ymin = min(ygrid(vertices_square(children,direction1)),ygrid(vertices_square(children,direction2)));
                    ymax = max(ygrid(vertices_square(children,direction1)),ygrid(vertices_square(children,direction2)));
                    indices = find(and(xgrid==xgrid(index),and(ygrid>ymin,ygrid<ymax)));
                    squares_vertex(indices,direction1) = squares_aux(index1);
                    squares_vertex(indices,direction2) = squares_aux(index1);
end
            while childrenQ(children,quadtree)
                children = quadtree(children,direction2);
            end
            %pointless
            %squares_vertex(index,index1) = children;
            
            children = quadtree(neighbour,direction2);
            switch cardinal_direction
                case {'N','S'}
                    xmin = min(xgrid(vertices_square(children,direction1)),xgrid(vertices_square(children,direction2)));
                    xmax = max(xgrid(vertices_square(children,direction1)),xgrid(vertices_square(children,direction2)));
                    indices = find(and(ygrid==ygrid(index),and(xgrid>xmin,xgrid<xmax)));
                    squares_vertex(indices,direction1) = squares_aux(index2);
                    squares_vertex(indices,direction2) = squares_aux(index2);
                case {'W','E'}
                    ymin = min(ygrid(vertices_square(children,direction1)),ygrid(vertices_square(children,direction2)));
                    ymax = max(ygrid(vertices_square(children,direction1)),ygrid(vertices_square(children,direction2)));
                    indices = find(and(xgrid==xgrid(index),and(ygrid>ymin,ygrid<ymax)));
                    squares_vertex(indices,direction1) = squares_aux(index2);
                    squares_vertex(indices,direction2) = squares_aux(index2);
            end
            while childrenQ(children,quadtree)
                children = quadtree(children,direction1);
            end
            %pointless
            %squares_vertex(index,index2) = children;
        else
            numTotGrid = numTotGrid+1;
            index = numTotGrid;
            vertices_index = vertices_square(current_square,:);
            xgrid(index) = mean(xgrid(vertices_index));
            ygrid(index) = mean(ygrid(vertices_index));
            l = side_length(current_square,vertices_square,xgrid);
            switch cardinal_direction
                case 'N'
                    ygrid(index) = ygrid(index) + l/2;
                case 'E'
                    xgrid(index) = xgrid(index) + l/2;
                case 'S'
                    ygrid(index) = ygrid(index) - l/2;
                case 'W'
                    xgrid(index) = xgrid(index) - l/2;
            end
            squares_vertex(index,index1) = neighbour;
            squares_vertex(index,index2) = neighbour;
        end
    else
        numTotGrid = numTotGrid+1;
        index = numTotGrid;
        vertices_index = vertices_square(current_square,:);
        xgrid(index) = mean(xgrid(vertices_index));
        ygrid(index) = mean(ygrid(vertices_index));
        l = side_length(current_square,vertices_square,xgrid);
        switch cardinal_direction
            case 'N'
                ygrid(index) = ygrid(index) + l/2;
            case 'E'
                xgrid(index) = xgrid(index) + l/2;
            case 'S'
                ygrid(index) = ygrid(index) - l/2;
            case 'W'
                xgrid(index) = xgrid(index) - l/2;
        end
    end
    switch cardinal_direction
        case 'N'
            index_N = index;
        case 'E'
            index_E = index;
        case 'S'
            index_S = index;
        case 'W'
            index_W = index;
    end
end
numTotGrid = numTotGrid+1;
index_C = numTotGrid;
vertices_index = vertices_square(current_square,:);
xgrid(index_C) = mean(xgrid(vertices_index));
ygrid(index_C) = mean(ygrid(vertices_index));

squares_vertex(index_C,1) = square_NW;
squares_vertex(index_C,2) = square_NE;
squares_vertex(index_C,3) = square_SW;
squares_vertex(index_C,4) = square_SE;


squares_vertex(index_N,3) = square_NW;
squares_vertex(index_N,4) = square_NE;

squares_vertex(index_E,1) = square_NE;
squares_vertex(index_E,3) = square_SE;


squares_vertex(index_S,1) = square_SW;
squares_vertex(index_S,2) = square_SE;

squares_vertex(index_W,2) = square_NW;
squares_vertex(index_W,4) = square_SW;

vertices_square(square_NW,1) = index_NW;
vertices_square(square_NW,2) = index_N;
vertices_square(square_NW,3) = index_C;
vertices_square(square_NW,4) = index_W;

vertices_square(square_NE,1) = index_N;
vertices_square(square_NE,2) = index_NE;
vertices_square(square_NE,3) = index_E;
vertices_square(square_NE,4) = index_C;

vertices_square(square_SW,1) = index_W;
vertices_square(square_SW,2) = index_C;
vertices_square(square_SW,3) = index_S;
vertices_square(square_SW,4) = index_SW;

vertices_square(square_SE,1) = index_C;
vertices_square(square_SE,2) = index_E;
vertices_square(square_SE,3) = index_SE;
vertices_square(square_SE,4) = index_S;
end