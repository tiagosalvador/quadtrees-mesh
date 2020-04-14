function iIntN = findneighbours(r,numDir,nu_list,numBdy,numInt,numTotGrid,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square)

% First quadrant (label 1)
% Fourth quadrant (label 2).
% Second quadrant (label 3).
% Third quadrant (label 4).


iIntN = zeros(numInt,4*numDir);
for k = 1:numDir
    node = (1:numTotGrid)';
    x0 = xgrid(node);
    y0 = ygrid(node);
    indices_inside_domain = find(insidedomain(x0,y0)>0);
    node = node(indices_inside_domain);
    
    
    nu_active = nu_list(:,k);
    lenNu = repmat(sqrt(nu_active(1,:).^2+nu_active(2,:).^2),2,1);
    nu_active = nu_active./lenNu;
    m = nu_active(2)/nu_active(1);
    
    neighbours_matrix = findneighboursaux(node,indices_inside_domain,m,nu_active,r,numBdy,numInt,x,y,xgrid,ygrid,xgrid2x,squares_boundary,insidedomain,squares_vertex,quadtree,vertices_square);
    
    iIntN(:,4*(k-1)+(1:4)) = neighbours_matrix((numBdy+1):(numBdy+numInt),:);
end
end
