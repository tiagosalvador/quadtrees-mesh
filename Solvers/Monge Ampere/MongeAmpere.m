function [F,JF] = MongeAmpere(u,Mbc,bc,f,numInt,iInt,iIntN,nu,numDir,numTot,numFr,x,y,eps)

% Compute all second derivatives.
d2u = SecondDeriv(u,x,y,iInt,iIntN,nu,numDir,numInt);

% MA operator
MA = max(d2u(:,1:numFr),eps) .* max(d2u(:,numFr+1:end),eps) + ...
    min(d2u(:,1:numFr),eps) + min(d2u(:,numFr+1:end),eps);

% Find the frame that gives the minimal value.
[MA,ac] = min(MA,[],2);

% Operator.
F = [Mbc*u-bc;-MA+f];

% The finite difference matrices corresponding to these active directions.
D2a = SecondDerivMatrix(x,y,iInt,iIntN,nu,numDir,numInt,numTot,ac);
D2b = SecondDerivMatrix(x,y,iInt,iIntN,nu,numDir,numInt,numTot,ac+numFr);

% Active second derivatives.
d2ua = d2u(sub2ind([numInt,2*numFr],(1:numInt)',ac));
d2ub = d2u(sub2ind([numInt,2*numFr],(1:numInt)',ac+numFr));

% Which part of operator is active.
con1 = double(d2ua>eps);
con2 = double(d2ub>eps);

% coefficients of second derivative matrices.
Ma = spdiags(con1.*max(d2ub,eps)+(1-con1),0,numInt,numInt);
Mb = spdiags(con2.*max(d2ua,eps)+(1-con2),0,numInt,numInt);

% Jacobian
JF = [Mbc;[-Ma*D2a - Mb*D2b]];

