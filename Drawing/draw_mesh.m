function draw_mesh(a,b,grid_parameters,mesh_parameters,domain_function)
figure;
indicator_boundary = @(x,y) domain_function(x,y);
insidedomain = @(x,y) domain_function(x,y)<0;
draw_quadtree(grid_parameters,mesh_parameters,insidedomain);
[q1,q2] = meshgrid(linspace(-2,2,1000));
hold on
contour(q1,q2,indicator_boundary(q1,q2),[0 0],'LineColor','black','LineWidth',2)
x = grid_parameters.x;
y = grid_parameters.y;
xB = grid_parameters.xB;
yB = grid_parameters.yB;
plot(x,y,'ok',xB,yB,'ob')
axis([a b a b])
hold off

end