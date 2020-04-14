function d2u = SecondDeriv(u,x,y,iInt,iN,nu,numDir,numInt)
%
% This function estimates the second directional derivative of a function u
% in the direction nu.  The approximation is monotone and relies on four
% neighbouring points that are nearly aligned with the directions \pm\nu.
%
% INPUTS:  u      A column vector containing the solution values at all
%                 discretisation points.
%          x,y    Column vectors containing the x and y coordinates of all
%                 discretisation points.  The boundary points should be
%                 listed first.
%          iInt   A column vector containing the indices of all interior
%                 points.
%          iN     A matrix of size length(iInt)X(4*numDir) containing the
%                 indices of the neighbours to use.  Each length(iInt)X4
%                 block corresponds to one direction.
%          nu     A 2XnumDir matrix containing the unit directions to use.
%          numDir The number of directions to consider.
%          numInt The number of interior points.
%
% OUTPUTS: d2u   A matrix with numDir columns containing the approximations
%                of the second derivatives at all of the interior points.
%

% Make the interior indices vector the same size as the neighbours vector.
iInt = repmat(iInt,1,4*numDir);

% Changes in (x,y) between neighbour and current point (column vectors).
dx = x(iN)-x(iInt);
dy = y(iN)-y(iInt);
if numInt==1
    dx = dx';
    dy = dy';
end;

% Augment.
nu1 = repmat(nu(1,:),4,1);
nu1 = repmat(nu1(:)',numInt,1);
nu2 = repmat(nu(2,:),4,1);
nu2 = repmat(nu2(:)',numInt,1);

% Sines and cosines made with nu, scaled by r.
C = nu1.*dx + nu2.*dy;
S = -nu2.*dx + nu1.*dy;

% Indices to find each neighbour.
i1 = 1:4:4*numDir-3;
i2 = 2:4:4*numDir-2;
i3 = 3:4:4*numDir-1;
i4 = 4:4:4*numDir;

% The approximation.
den = ((C(:,i3).*S(:,i4)-C(:,i4).*S(:,i3)).*(C(:,i1).^2.*S(:,i2)-C(:,i2).^2.*S(:,i1))-...
    (C(:,i1).*S(:,i2)-C(:,i2).*S(:,i1)).*(C(:,i3).^2.*S(:,i4)-C(:,i4).^2.*S(:,i3)));
co1 = 2*S(:,i2).*(C(:,i3).*S(:,i4)-C(:,i4).*S(:,i3))./den;
co2 = -2*S(:,i1).*(C(:,i3).*S(:,i4)-C(:,i4).*S(:,i3))./den;
co3 = -2*S(:,i4).*(C(:,i1).*S(:,i2)-C(:,i2).*S(:,i1))./den;
co4 = 2*S(:,i3).*(C(:,i1).*S(:,i2)-C(:,i2).*S(:,i1))./den;
co0 = -(co1+co2+co3+co4);
u1 = u(iN(:,i1));
u2 = u(iN(:,i2));
u3 = u(iN(:,i3));
u4 = u(iN(:,i4));
u0 = u(iInt(:,i1));
if numInt==1
    u1 = u1'; u2 = u2'; u3 = u3'; u4 = u4'; u0 = u0';
end;
d2u = co1.*u1 + co2.*u2 + co3.*u3 + co4.*u4 + co0.*u0;
