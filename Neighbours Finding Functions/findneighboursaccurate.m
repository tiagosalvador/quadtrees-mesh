function [ind, axx, ayy, axy] = findneighboursaccurate(numBdy,numTot,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square,indicator_boundary)

% we assume that the boundary only intersects an edge of any square
% at most once

ind = zeros(numTot,9);
axx = zeros(numTot,9);
ayy = zeros(numTot,9);
axy = zeros(numTot,9);
% The ind is structured in the following way: each row
% represents the following information of a interior point:
%       1st column - N grid neighbour
%       2nd column - E grid neighbour
%       3rd column - S grid neighbour
%       4th column - W grid neighbour
%       5th column - NW grid neighbour
%       6th column - NE grid neighbour
%       7th column - SE grid neighbour
%       8th column - SW grid neighbour
%       9th column - extra grid neighbour

for node = 1:numTotGrid
        x0 = xgrid(node);
        y0 = ygrid(node);
        if insidedomain(x0,y0)
            squares = squares_vertex(node,1:4);
            depths = depth_level(squares,quadtree);
            while length(unique(depths)) > 1
                indices = depths > min(depths);
                squares(indices) = quadtree(squares(indices),5);
                depths(indices) = depth_level(squares(indices),quadtree);
            end
            vertices_index = unique(vertices_square(squares,:));
            
            if length(unique(squares)) == 4
                % N grid neighbour
                current_square = squares(1);
                index = vertices_square(current_square,2);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),1) = xgrid2x(index);
                else
                    ind(xgrid2x(node),1) = squares_boundary(current_square,2);
                end
                % E grid neighbour
                current_square = squares(2);
                index = vertices_square(current_square,3);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),2) = xgrid2x(index);
                else
                    ind(xgrid2x(node),2) = squares_boundary(current_square,3);
                end
                % S grid neighbour
                current_square = squares(4);
                index = vertices_square(current_square,4);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),3) = xgrid2x(index);
                else
                    ind(xgrid2x(node),3) = squares_boundary(current_square,4);
                end
                % W grid neighbour
                current_square = squares(3);
                index = vertices_square(current_square,1);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),4) = xgrid2x(index);
                else
                    ind(xgrid2x(node),4) = squares_boundary(current_square,1);
                end
                % NW grid
                current_square = squares(1);
                index = vertices_square(current_square,1);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),5) = xgrid2x(index);
                end
                % NE grid
                current_square = squares(2);
                index = vertices_square(current_square,2);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),6) = xgrid2x(index);
                end
                % SE grid
                current_square = squares(4);
                index = vertices_square(current_square,3);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),7) = xgrid2x(index);
                end
                % SW grid
                current_square = squares(3);
                index = vertices_square(current_square,4);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),8) = xgrid2x(index);
                end
            elseif length(unique(squares)) == 2
                
                if squares(1) == squares(2)
                    if childrenQ(squares(1),quadtree)
                        situation = 1;
                    else
                        situation = 2;
                    end
                else
                    if childrenQ(squares(1),quadtree)
                        situation  = 3;
                    else
                        situation = 4;
                    end
                end
                
                % N grid neighbour
                if not(situation == 2)
                    if situation == 1
                        current_square = quadtree(squares(1),3);
                    else
                        current_square = squares(1);
                    end
                    index = vertices_square(current_square,2);
                    if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                        ind(xgrid2x(node),1) = xgrid2x(index);
                    else
                        ind(xgrid2x(node),1) = squares_boundary(current_square,2);
                    end
                end
                % E grid neighbour
                if not(situation == 3)
                    if situation == 4
                        current_square = quadtree(squares(2),1);
                    else
                        current_square = squares(1);
                    end
                    index = vertices_square(current_square,3);
                    if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                        ind(xgrid2x(node),2) = xgrid2x(index);
                    else
                        ind(xgrid2x(node),2) = squares_boundary(current_square,3);
                    end
                end
                % S grid neighbour
                if not(situation == 1)
                    if situation == 2
                        current_square = quadtree(squares(3),2);
                    else
                        current_square = squares(2);
                    end
                    index = vertices_square(current_square,4);
                    if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                        ind(xgrid2x(node),3) = xgrid2x(index);
                    else
                        ind(xgrid2x(node),3) = squares_boundary(current_square,4);
                    end
                end
                % W grid neighbour
                if not(situation == 4)
                    if situation == 3
                        current_square = quadtree(squares(1),4);
                    else
                        current_square = squares(3);
                    end
                    index = vertices_square(current_square,1);
                    if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                        ind(xgrid2x(node),4) = xgrid2x(index);
                    else
                        ind(xgrid2x(node),4) = squares_boundary(current_square,1);
                    end
                end
                % NW grid
                current_square = squares(1);
                index = vertices_square(current_square,1);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),5) = xgrid2x(index);
                end
                % NE grid
                current_square = squares(2);
                index = vertices_square(current_square,2);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),6) = xgrid2x(index);
                end
                % SE grid
                current_square = squares(4);
                index = vertices_square(current_square,3);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),7) = xgrid2x(index);
                end
                % SW grid
                current_square = squares(3);
                index = vertices_square(current_square,4);
                if indicator_boundary(xgrid(index),ygrid(index)) <= 0
                    ind(xgrid2x(node),8) = xgrid2x(index);
                end
                
                
            else
                disp('this should not happen')
            end
            
            if not(prod(indicator_boundary(xgrid(vertices_index),ygrid(vertices_index)) <= 0))
                % easy thing: pick whatever nodes from the boundary
                % with no criteria. This can be improved
                
                %                 indices_left = find(ind(xgrid2x(node),:)==0);
                %                 indices_boundary = setdiff(unique(squares_boundary(squares,5:end)),0);
                %                 ind(xgrid2x(node),indices_left) = indices_boundary(1:length(indices_left));
                %
                %                 aux = ind(xgrid2x(node),:);
                %
                %
                %                 M = [(x(aux)'-x0); (y(aux)'-y0); (x(aux)'-x0).^2/2; (y(aux)'-y0).^2/2;...
                %                     (x(aux)'-x0).*(y(aux)'-y0); (x(aux)'-x0).^3/6; (y(aux)'-y0).^3/6;...
                %                     (x(aux)'-x0).^2.*(y(aux)'-y0)/2; (x(aux)'-x0).*(y(aux)'-y0).^2/2];
                %
                %                 axx(xgrid2x(node),:) = M\[zeros(2,1);1;zeros(6,1)];
                %                 ayy(xgrid2x(node),:) = M\[zeros(3,1);1;zeros(5,1)];
                %                 axy(xgrid2x(node),:) = M\[zeros(4,1);1;zeros(4,1)];
                %
                xmin = min(min(xgrid(vertices_square(squares,:))));
                xmax = max(max(xgrid(vertices_square(squares,:))));
                ymin = min(min(ygrid(vertices_square(squares,:))));
                ymax = max(max(ygrid(vertices_square(squares,:))));
                aux = find(and(x>=xmin,and(x<=xmax,and(y>=ymin,y<=ymax))));
                
                ind(xgrid2x(node),1:length(aux)) = aux;
                
                
                M = [(x(aux)'-x0); (y(aux)'-y0); (x(aux)'-x0).^2/2; (y(aux)'-y0).^2/2;...
                    (x(aux)'-x0).*(y(aux)'-y0); (x(aux)'-x0).^3/6; (y(aux)'-y0).^3/6;...
                    (x(aux)'-x0).^2.*(y(aux)'-y0)/2; (x(aux)'-x0).*(y(aux)'-y0).^2/2];
                %                 axx(xgrid2x(node),1:length(aux)) = leastnormsolution(M,[zeros(2,1);1;zeros(6,1)]);
                %                 ayy(xgrid2x(node),1:length(aux)) = leastnormsolution(M,[zeros(3,1);1;zeros(5,1)]);
                %                 axy(xgrid2x(node),1:length(aux)) = leastnormsolution(M,[zeros(4,1);1;zeros(4,1)]);
                
                if rank(M*M') < 9
                      squaresN = north_neighbour_vector(squares',quadtree);
                      squaresW = west_neighbour_vector(squares',quadtree);
                      squaresS = south_neighbour_vector(squares',quadtree);
                      squaresE = east_neighbour_vector(squares',quadtree);
                      squares = union(squaresN,union(squaresW,union(squaresS,union(squaresE,squares))));
                      squares = unique(squares);
                      squares = squares(squares>0);
                      xmin = min(min(xgrid(vertices_square(squares,:))));
                      xmax = max(max(xgrid(vertices_square(squares,:))));
                      ymin = min(min(ygrid(vertices_square(squares,:))));
                      ymax = max(max(ygrid(vertices_square(squares,:))));
                      aux = find(and(x>=xmin,and(x<=xmax,and(y>=ymin,y<=ymax))));
                      
                      ind(xgrid2x(node),1:length(aux)) = aux;
                      M = [(x(aux)'-x0); (y(aux)'-y0); (x(aux)'-x0).^2/2; (y(aux)'-y0).^2/2;...
                          (x(aux)'-x0).*(y(aux)'-y0); (x(aux)'-x0).^3/6; (y(aux)'-y0).^3/6;...
                          (x(aux)'-x0).^2.*(y(aux)'-y0)/2; (x(aux)'-x0).*(y(aux)'-y0).^2/2];
                     if rank(M*M')<9
                         disp('oh no')
                     end
                end
                
                Maux =  M'*inv(M*M');
                axx(xgrid2x(node),1:length(aux)) = Maux*[zeros(2,1);1;zeros(6,1)];
                ayy(xgrid2x(node),1:length(aux)) = Maux*[zeros(3,1);1;zeros(5,1)];
                axy(xgrid2x(node),1:length(aux)) = Maux*[zeros(4,1);1;zeros(4,1)];
                
%                 axx(xgrid2x(node),1:length(aux)) = M\[zeros(2,1);1;zeros(6,1)];
%                 ayy(xgrid2x(node),1:length(aux)) = M\[zeros(3,1);1;zeros(5,1)];
%                 axy(xgrid2x(node),1:length(aux)) = M\[zeros(4,1);1;zeros(4,1)];
            else
                dx = side_length(squares(1),vertices_square,xgrid);
                if length(unique(squares)) == 4
                    axx(xgrid2x(node),[2 4]) = 1/dx^2;
                    ayy(xgrid2x(node),[1 3]) = 1/dx^2;
                    axy(xgrid2x(node),[5 7]) = -1/(4*dx^2);
                    axy(xgrid2x(node),[6 8]) = 1/(4*dx^2);
                elseif length(unique(squares)) == 2
                    if situation <= 2
                        axx(xgrid2x(node),[2 4]) = 4/dx^2;
                        ayy(xgrid2x(node),[2 4]) = -1/dx^2;
                        ayy(xgrid2x(node),4:8) = 1/(2*dx^2);
                    else
                        axx(xgrid2x(node),[1 3]) = -1/dx^2;
                        axx(xgrid2x(node),4:8) = 1/(2*dx^2);
                        ayy(xgrid2x(node),[1 3]) = 4/dx^2;
                    end
                    axy(xgrid2x(node),[5 7]) = -1/(2*dx^2);
                    axy(xgrid2x(node),[6 8]) = 1/(2*dx^2);
                end
            end
            
        end
end

ind(1:numBdy,:) = [];
axx(1:numBdy,:) = [];
ayy(1:numBdy,:) = [];
axy(1:numBdy,:) = [];