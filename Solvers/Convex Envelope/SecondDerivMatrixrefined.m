function D2 = SecondDerivMatrixrefined(x,y,iInt,iN,nu1,nu2,numDir,numInt,numTot,iNu)
%
% This function provides the finite difference matrix for the second
% directional derivatives in given active directions.
%
% INPUTS:  x,y    Column vectors containing the x and y coordinates of all
%                 discretisation points.  The boundary points should be
%                 listed first.
%          iInt   A column vector containing the indices of all interior
%                 points.
%          iN     A matrix of size numIntX(4*numDir) containing the
%                 indices of the neighbours to use.  Each numIntX4
%                 block corresponds to one direction.
%          nu     A 2XnumDir matrix containing the unit directions to use.
%          numDir The number of directions to consider.
%          numInt The number of interior points.
%          numTot The total number of discretisation points.
%          iNu    A column vector of length numInt giving the index of the
%                 direction that is active at each interior point.  All
%                 entries should be integers between one and numDir.
%          
% OUTPUTS: D2     A sparse matrix of size numInt X numTot that is the
%                 finite difference matrix for the second directional
%                 derivatives in the active directions.
%

% Choose the appropriate neighbours.
iUse1 = sub2ind([numInt,4*numDir],(1:numInt)',4*(iNu-1)+1);
iUse2 = sub2ind([numInt,4*numDir],(1:numInt)',4*(iNu-1)+2);
iUse3 = sub2ind([numInt,4*numDir],(1:numInt)',4*(iNu-1)+3);
iUse4 = sub2ind([numInt,4*numDir],(1:numInt)',4*iNu);
iN = [iN(iUse1),iN(iUse2),iN(iUse3),iN(iUse4)];

% Make the interior indices vector the same size as the neighbours vector.
iInt = repmat(iInt,1,4);

% Changes in (x,y) between neighbour and current point (column vectors).
dx = x(iN)-x(iInt);
dy = y(iN)-y(iInt);

% Augment.
nu1 = repmat(nu1',1,4);
nu2 = repmat(nu2',1,4);

% Sines and cosines made with nu, scaled by r.
C = nu1.*dx + nu2.*dy;
S = -nu2.*dx + nu1.*dy;

% The approximation.
den = ((C(:,3).*S(:,4)-C(:,4).*S(:,3)).*(C(:,1).^2.*S(:,2)-C(:,2).^2.*S(:,1))-...
    (C(:,1).*S(:,2)-C(:,2).*S(:,1)).*(C(:,3).^2.*S(:,4)-C(:,4).^2.*S(:,3)));
co1 = 2*S(:,2).*(C(:,3).*S(:,4)-C(:,4).*S(:,3))./den;
co2 = -2*S(:,1).*(C(:,3).*S(:,4)-C(:,4).*S(:,3))./den;
co3 = -2*S(:,4).*(C(:,1).*S(:,2)-C(:,2).*S(:,1))./den;
co4 = 2*S(:,3).*(C(:,1).*S(:,2)-C(:,2).*S(:,1))./den;
co0 = -(co1+co2+co3+co4);

% The second derivative matrix.
D2 = sparse(1:numInt,iN(:,1),co1,numInt,numTot) + ...
    sparse(1:numInt,iN(:,2),co2,numInt,numTot) + ...
    sparse(1:numInt,iN(:,3),co3,numInt,numTot) + ...
    sparse(1:numInt,iN(:,4),co4,numInt,numTot) + ...
    sparse(1:numInt,iInt(:,1),co0,numInt,numTot);

