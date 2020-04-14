function r = gConvEnv(x,y)
phi = pi/6;
Rx = cos(-phi)*x-sin(-phi)*y;
Ry = sin(-phi)*x + cos(-phi)*y;
r = min(min(sqrt((Rx+0.5).^2+Ry.^2),...
    sqrt((Rx-0.5).^2+Ry.^2)),0.5);
end

