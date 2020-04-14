function r = domain_conv_env(x,y)
phi = pi/6;
Rx = cos(-phi)*x-sin(-phi)*y;
Ry = sin(-phi)*x + cos(-phi)*y;
g1 = (Rx+0.5).^2+Ry.^2-0.5^2;
g2 = (Rx-0.5).^2+Ry.^2-0.5^2;
g3 = max(abs(Rx),abs(Ry))-0.5;
r = min(g1,min(g2,g3));
end