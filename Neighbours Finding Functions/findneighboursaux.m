function neighbours_matrix = findneighboursaux(nodes_to_update,indices_inside_domain,m,nu_active,r,numBdy,numInt,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square)
neighbours_matrix = zeros(numBdy+numInt,4);    
node = nodes_to_update;
x0 = xgrid(node);
y0 = ygrid(node);

direction = 1;
while direction<= 2
    if direction == 1
        if m > 0
            current_square = squares_vertex(indices_inside_domain,2);
            edge1 = 'E';
            edge2 = 'N';
            vertice_label = 2;
        else
            current_square = squares_vertex(indices_inside_domain,4);
            edge1 = 'S';
            edge2 = 'E';
            vertice_label = 3;
        end
    else
        node = nodes_to_update;
        x0 = xgrid(node);
        y0 = ygrid(node);
        if m > 0
            current_square = squares_vertex(indices_inside_domain,3);
            edge1 = 'W';
            edge2 = 'S';
            vertice_label = 4;
        else
            current_square = squares_vertex(indices_inside_domain,1);
            edge1 = 'N';
            edge2 = 'W';
            vertice_label = 1;
        end
    end
    %indices_inside_domain = indices_inside_domain(current_square > 0);
    %current_square = current_square(current_square > 0);
    go = ~isempty(current_square);
    while go
        go = 0;
        
        x1 = xgrid(vertices_square(current_square,vertice_label));
        y1 = ygrid(vertices_square(current_square,vertice_label));
        
        
        switch edge1
            case {'E','W'}
                x_aux = x1;
                y_aux = m*(x1-x0)+y0;
            case {'N','S'}
                x_aux = (y1-y0)/m+x0;
                y_aux = y1;
        end
        
        switch edge1
            case 'E'
                check = and(ygrid(vertices_square(current_square,3)) <= y_aux,...
                    y_aux <= ygrid(vertices_square(current_square,2)));
            case 'S'
                check = and(xgrid(vertices_square(current_square,4)) <= x_aux,...
                    x_aux <= xgrid(vertices_square(current_square,3)));
            case 'W'
                check = and(ygrid(vertices_square(current_square,4)) <= y_aux,...
                    y_aux <= ygrid(vertices_square(current_square,1)));
            case 'N'
                check = and(xgrid(vertices_square(current_square,1)) <= x_aux,...
                    x_aux <= xgrid(vertices_square(current_square,2)));
        end
        
        future_square = zeros(size(current_square));
        aux = zeros(length(current_square),2);
        
        switch edge1
            case 'E'
                aux(check,:) = repmat([2 3],nnz(check),1);
                future_square(check) = east_neighbour_vector(current_square(check),quadtree,vertices_square,x_aux(check),y_aux(check),xgrid,ygrid);
                aux(and(check,future_square>0),:) = repmat([1 4],nnz(and(check,future_square>0)),1);
            case 'S'
                aux(check,:) = repmat([4 3],nnz(check),1);
                future_square(check) = south_neighbour_vector(current_square(check),quadtree,vertices_square,x_aux(check),y_aux(check),xgrid,ygrid);
                aux(and(check,future_square>0),:) = repmat([1 2],nnz(and(check,future_square>0)),1);
            case 'W'
                aux(check,:) = repmat([1 4],nnz(check),1);
                future_square(check) = west_neighbour_vector(current_square(check),quadtree,vertices_square,x_aux(check),y_aux(check),xgrid,ygrid);
                aux(and(check,future_square>0),:) = repmat([2 3],nnz(and(check,future_square>0)),1);
            case 'N'
                aux(check,:) = repmat([1 2],nnz(check),1);
                future_square(check) = north_neighbour_vector(current_square(check),quadtree,vertices_square,x_aux(check),y_aux(check),xgrid,ygrid);
                aux(and(check,future_square>0),:) = repmat([4 3],nnz(and(check,future_square>0)),1);
        end
        
        switch edge2
            case {'E','W'}
                x_aux(not(check)) = x1(not(check));
                y_aux(not(check)) = m*(x1(not(check))-x0(not(check)))+y0(not(check));
            case {'N','S'}
                x_aux(not(check)) = (y1(not(check))-y0(not(check)))/m+x0(not(check));
                y_aux(not(check)) = y1(not(check));
        end
        switch edge2
            case 'E'
                aux(not(check),:) = repmat([2 3],nnz(not(check)),1);
                future_square(not(check)) = east_neighbour_vector(current_square(not(check)),quadtree,vertices_square,x_aux(not(check)),y_aux(not(check)),xgrid,ygrid);
                aux(and(not(check),future_square>0),:) = repmat([1 4],nnz(and(not(check),future_square>0)),1);
            case 'S'
                aux(not(check),:) = repmat([4 3],nnz(not(check)),1);
                future_square(not(check)) = south_neighbour_vector(current_square(not(check)),quadtree,vertices_square,x_aux(not(check)),y_aux(not(check)),xgrid,ygrid);
                aux(and(not(check),future_square>0),:) = repmat([1 2],nnz(and(not(check),future_square>0)),1);
            case 'W'
                aux(not(check),:) = repmat([1 4],nnz(not(check)),1);
                future_square(not(check)) = west_neighbour_vector(current_square(not(check)),quadtree,vertices_square,x_aux(not(check)),y_aux(not(check)),xgrid,ygrid);
                aux(and(not(check),future_square>0),:) = repmat([2 3],nnz(and(not(check),future_square>0)),1);
            case 'N'
                aux(not(check),:) = repmat([1 2],nnz(not(check)),1);
                future_square(not(check)) = north_neighbour_vector(current_square(not(check)),quadtree,vertices_square,x_aux(not(check)),y_aux(not(check)),xgrid,ygrid);
                aux(and(not(check),future_square>0),:) = repmat([4 3],nnz(and(not(check),future_square>0)),1);
        end
        active_square = zeros(size(future_square));
        active_square(future_square > 0) = future_square(future_square > 0);
        active_square(not(future_square > 0)) = current_square(not(future_square > 0));
        % previously it was vectorized, but apparently that is slower
        update_these = sub2ind(size(vertices_square),active_square,aux(:,1));
        dist1 = sqrt((xgrid(vertices_square(update_these))-xgrid(node)).^2+(ygrid(vertices_square(update_these))-ygrid(node)).^2);
        update_these = sub2ind(size(vertices_square),active_square,aux(:,2));
        dist2 = sqrt((xgrid(vertices_square(update_these))-xgrid(node)).^2+(ygrid(vertices_square(update_these))-ygrid(node)).^2);
        
        for kk = 1:2
            index = aux(:,kk);
            update_these = sub2ind(size(vertices_square),active_square,index);
            i = vertices_square(update_these);
            %i = diag(vertices_square(active_square,index));
            rr = sqrt((xgrid(i)-xgrid(node)).^2+(ygrid(i)-ygrid(node)).^2);
            C = (nu_active(1)*(xgrid(i)-xgrid(node)) + nu_active(2)*(ygrid(i)-ygrid(node))) ./ rr;
            S = (-nu_active(2)*(xgrid(i)-xgrid(node)) + nu_active(1)*(ygrid(i)-ygrid(node))) ./ rr;
            label = zeros(size(S));
            label(and(C>0,S>=0)) = 1;
            label(and(C>=0,S<0)) = 2;
            label(and(C<=0,S>0)) = 3;
            label(and(C<0,S<=0)) = 4;
            check = and((rr <= r),insidedomain(xgrid(i),ygrid(i)));
            i = xgrid2x(i(check));
            %icurrent = diag(neighbours_matrix(xgrid2x(node(check)),label(check)))';
            update_these = sub2ind(size(neighbours_matrix),xgrid2x(node(check)),label(check)');
            icurrent = neighbours_matrix(update_these);
            neighbours_matrix(update_these) = update_neighbours_matrix_vector(xgrid2x(node(check)),i,icurrent,nu_active,x,y);
        end
        
        neighbours_matrix = update_boundary_vector(node,neighbours_matrix,current_square,squares_boundary,nu_active,x,y,xgrid2x,xgrid,ygrid,r);
        x0 = x_aux;
        y0 = y_aux;
        
        dist = sqrt((x0-xgrid(node)).^2+(y0-ygrid(node)).^2);
        %check = not(or(and(r<dist1,r<dist2),not(insidedomain(x0,y0))));
        check = and(insidedomain(x0,y0),dist<r);
        node = node(check);
        current_square = future_square(check);
        x0 = x0(check);
        y0 = y0(check);
        
        check = current_square>0;
        x0 = x0(check);
        y0 = y0(check);
        node = node(check);
        current_square = current_square(check);
        go = ~isempty(current_square);
    end
    direction = direction + 1;
end
end