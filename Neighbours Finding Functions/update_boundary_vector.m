function neighbours_matrix = update_boundary_vector(node_vector,neighbours_matrix,current_square_vector,square_boundary,nu,x,y,xgrid2x,xgrid,ygrid,r)

for kk = 1:size(square_boundary,2)
    
    boundary = square_boundary(current_square_vector,kk);
    node = node_vector(boundary>0);
    boundary = boundary(boundary>0);
    rr = sqrt((x(boundary)-xgrid(node)).^2+(y(boundary)-ygrid(node)).^2);
    node = node(rr <= r);
    boundary = boundary(rr <= r);
    rr = rr(rr <= r);
    if and(kk>4,isempty(boundary))
        break
    else
        C = (nu(1)*(x(boundary)-xgrid(node)) + nu(2)*(y(boundary)-ygrid(node))) ./ rr;
        S = (-nu(2)*(x(boundary)-xgrid(node)) + nu(1)*(y(boundary)-ygrid(node))) ./ rr;
        label = zeros(size(S));
        label(and(C>0,S>=0)) = 1;
        label(and(C>=0,S<0)) = 2;
        label(and(C<=0,S>0)) = 3;
        label(and(C<0,S<=0)) = 4;
        icurrent = diag(neighbours_matrix(xgrid2x(node),label))';
        update_these = sub2ind(size(neighbours_matrix),xgrid2x(node),label');
        neighbours_matrix(update_these) = update_neighbours_matrix_vector(xgrid2x(node),boundary,icurrent,nu,x,y);
    end
end


% for kk=1:length(node_vector)
%     node = node_vector(kk);
%     current_square = current_square_vector(kk);
%     if current_square
%         for i = 1:length(square_boundary(current_square,:))
%             boundary = square_boundary(current_square,i);
%             if boundary
%                 rr = sqrt((x(boundary)-xgrid(node)).^2+(y(boundary)-ygrid(node)).^2);
%                 if rr <= r
%                     C = (nu(1)*(x(boundary)-xgrid(node)) + nu(2)*(y(boundary)-ygrid(node))) ./ rr;
%                     S = (-nu(2)*(x(boundary)-xgrid(node)) + nu(1)*(y(boundary)-ygrid(node))) ./ rr;                 
%                     if and(C>0,S>=0)
%                         label = 1;
%                     elseif and(C>=0,S<0)
%                         label = 2;
%                     elseif and(C<=0,S>0)
%                         label = 3;
%                     elseif and(C<0,S<=0)
%                         label = 4;
%                     end
%                     icurrent = neighbours_matrix(xgrid2x(node),label);
%                     neighbours_matrix(xgrid2x(node),label) = update_neighbours_matrix(xgrid2x(node),boundary,icurrent,nu,x,y,0);
%                 end
%             end
%         end
%     end
% end
end
