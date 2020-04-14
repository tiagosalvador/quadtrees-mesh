function [F,JF] = MongeAmpereFiltered(u,Mbc,bc,f,numInt,iInt,iIntN,nu,numDir,numTot,numFr,x,y,eps,Dxx,Dyy,Dxy,epsilon)


% Compute all second derivatives.
d2u = SecondDeriv(u,x,y,iInt,iIntN,nu,numDir,numInt);

% MA operator
MA = max(d2u(:,1:numFr),eps) .* max(d2u(:,numFr+1:end),eps) + ...
    min(d2u(:,1:numFr),eps) + min(d2u(:,numFr+1:end),eps);

% Find the frame that gives the minimal value.
[MA,ac] = min(MA,[],2);

% Values of standard derivatives.
uxx = Dxx*u;
uyy = Dyy*u;
uxy = Dxy*u;

% MA operator.
Fmon = -MA + f;
Fac = -(uxx.*uyy-uxy.^2) + f;

% Filter.
S = @(x)(abs(x)<=1).*x + (x>1).*max(2-x,0) + (x<-1).*min(-2-x,0);
dS = @(x)double(abs(x)<1) - double(and(abs(x)>1,abs(x)<2));

% Filtered operator.
Ffil = Fmon + epsilon * S((Fac-Fmon)/epsilon);

% Operator.
F = [Mbc*u-bc;Ffil];

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
Mxx = spdiags(uyy,0,numInt,numInt);
Myy = spdiags(uxx,0,numInt,numInt);
Mxy = spdiags(2*uxy,0,numInt,numInt);

% The derivative of the filter.
dS = dS((Fac-Fmon) / epsilon);
dS = max(dS,0);

% The weights attached to each gradient.
maskM = spdiags(1-dS,0,numInt,numInt);
maskA = spdiags(dS,0,numInt,numInt);

% Jacobian of MA operator.
Jac = maskM*(-Ma*D2a - Mb*D2b) + maskA*(-Mxx*Dxx-Myy*Dyy+Mxy*Dxy);

% Jacobian
JF = [Mbc;Jac];

