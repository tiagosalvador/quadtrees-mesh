function r = domain_union_circles(x,y)
% computation domain should be [-1.5,1.5]^2
circle1 = (x-0.25).^2+y.^2-1^2;
circle2 = (x+0.5).^2+y.^2-0.7^2;
r = min(circle1,circle2);
end