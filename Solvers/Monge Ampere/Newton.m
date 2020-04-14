function [x,res,f] = Newton(fun,x0,tol,alphaMin)

[f,df] = fun(x0);
% h = 1/sqrt(length(x0));
% res = norm(f(:),2)*h
res = norm(f(:),inf);
% alpha = 1;
alpha= 1;
x = x0;
its = 0;
while and(res>tol,alpha>alphaMin)
    its = its+1;
    alpha = 1;
    dX = df\f;
    [f,df] = fun(x-dX);
%     resNew = norm(f(:),2)*h;
    resNew = norm(f(:),inf);
    while and(resNew>res,alpha>alphaMin)
        alpha = 0.5*alpha;
        [f,df] = fun(x-alpha*dX);
%         resNew = norm(f(:),2)*h;
        resNew = norm(f(:),inf);
    end
    if alpha>alphaMin
        x = x-alpha*dX;
        res = resNew;
    end
end
res;
end