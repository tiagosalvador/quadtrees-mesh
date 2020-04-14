function r = domain_clover(x,y)
% computation domain should be [-1.5,1.5]^2
r = min((x-1/2).^2+5*(y-1/4).^2-1/2,min(5*(x+1/4).^2+(y-1/4).^2-1/2,min(5*(x-1/4).^2+(y+1/4).^2-1/2,(x+1/2).^2+5*(y+1/4).^2-1/2)));
end
