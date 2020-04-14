function [u,G] = CEActive(u,Mbc,bc,g,numInt,iInt,iBdy,iIntN,nu,numDir,numTot,x,y,tol)


% Compute all second derivatives.
d2u = SecondDeriv(u,x,y,iInt,iIntN,nu,numDir,numInt);

% Minimal eigenvalue
[lam1,ac] = min(d2u,[],2);

% Is the obstacle active?
obs = (u(iInt)-g>-lam1);

% Residual 
res1 = max(-lam1,u(iInt)-g);
res2 = u(iBdy)-bc;
res = norm([res1;res2],inf);

% tri = delaunay(x,y);

% Iterate.
while res > tol

% The finite difference matrices corresponding to these active directions.
D2 = SecondDerivMatrix(x,y,iInt,iIntN,nu,numDir,numInt,numTot,ac);

% Identity matrix.
Id = sparse(1:numInt,iInt,1,numInt,numTot);

% Coefficients of matrices.
coD2 = spdiags(1-obs,0,numInt,numInt);
coId = spdiags(obs,0,numInt,numInt);

% Matrix
M = [Mbc;coD2*D2+coId*Id];

% Rhs.
rhs = [zeros(numInt,1)];
rhs(obs) = g(obs);

% Solve.
u = M\[bc;rhs];

% Compute all second derivatives.
d2u = SecondDeriv(u,x,y,iInt,iIntN,nu,numDir,numInt);

% Minimal eigenvalue.
[lam1,ac] = min(d2u,[],2);

% Is the obstacle active?
obs = (u(iInt)-g>-lam1);

% trisurf(tri,x,y,u); drawnow;

% Residual
res1 = max(-lam1,u(iInt)-g);
res2 = u(iBdy)-bc;
res = norm([res1;res2],inf);

end;
G = [res2;min(abs(lam1),abs(u(iInt)-g))];
end
