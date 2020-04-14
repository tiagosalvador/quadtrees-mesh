function r = update_neighbours_matrix_vector(node,i,icurrent,nu,x,y)

r = icurrent;
r(icurrent == 0) = i(icurrent == 0);
icurrent(icurrent==0) = i(icurrent == 0);
if nnz(icurrent) > 0
    ri = sqrt((x(node)-x(i)).^2+(y(node)-y(i)).^2);
    Si = (-nu(2)*(x(i)-x(node)) + nu(1)*(y(i)-y(node))) ./ ri;
    ricurrent = sqrt((x(node)-x(icurrent)).^2+(y(node)-y(icurrent)).^2);
    Sicurrent = (-nu(2)*(x(icurrent)-x(node)) + nu(1)*(y(icurrent)-y(node))) ./ ricurrent;
    aux = find(or(abs(Si) < abs(Sicurrent),and(Si==Sicurrent,ri<ricurrent)));
    r(aux) = i(aux);
end
end