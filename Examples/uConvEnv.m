function r = uConvEnv(x,y)
phi = pi/6;
Rx = cos(-phi)*x-sin(-phi)*y;
Ry = sin(-phi)*x + cos(-phi)*y;
r = (abs(Rx)>0.5).*min(sqrt((Rx+0.5).^2+Ry.^2),...
    sqrt((Rx-0.5).^2+Ry.^2)) + ...
    (abs(Rx)<=0.5).*abs(Ry);
end

