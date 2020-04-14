function [u,G] = CEActiverefined(refinement_level,u,Mbc,bc,g,numInt,iInt,iBdy,iIntN,nu,numDir,numTot,x,y,tol)

d2u = NaN(numInt,max(numDir));
% Compute all second derivatives.
for kk=1:length(numDir)
    d2u(refinement_level==kk,1:numDir(kk)) = SecondDeriv(u,x,y,iInt(refinement_level==kk),iIntN(refinement_level==kk,1:4*numDir(kk)),nu{kk},numDir(kk),length(iInt(refinement_level==kk)));
end
% Minimal eigenvalue
[lam1,ac] = min(d2u,[],2);

% Is the obstacle active?
obs = (u(iInt)-g>-lam1);

% Residual 
res1 = max(-lam1,u(iInt)-g);
res2 = u(iBdy)-bc;
res = norm([res1;res2],inf);

% tri = delaunay(x,y);

nu1 = zeros(1,numInt);
nu2 = zeros(1,numInt);
% Iterate.
while res > tol

% The finite difference matrices corresponding to these active directions.
% Chosse the appropriate directions.

for kk=1:length(numDir)
    nu1(refinement_level==kk) = nu{kk}(1,ac(refinement_level==kk));
    nu2(refinement_level==kk) = nu{kk}(2,ac(refinement_level==kk));
end
D2 = SecondDerivMatrixrefined(x,y,iInt,iIntN,nu1,nu2,max(numDir),numInt,numTot,ac);

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
d2u = NaN(numInt,max(numDir));
% Compute all second derivatives.
for kk=1:length(numDir)
    d2u(refinement_level==kk,1:numDir(kk)) = SecondDeriv(u,x,y,iInt(refinement_level==kk),iIntN(refinement_level==kk,1:4*numDir(kk)),nu{kk},numDir(kk),length(iInt(refinement_level==kk)));
end

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
