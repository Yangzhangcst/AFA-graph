function W = OLRSCGRAPH_n(feature,para)
 
    epochs = para.epochs;
    d = para.d;
    rho = para.rho;

    Z = feature';
    [pp, nn] = size(Z);
    lambda1 = 1;

    lambda2_base = 1/sqrt(pp);
    MM = zeros(pp, d);
    AA = zeros(d, d);
    BB = zeros(pp, d);
    UU = zeros(nn, d);
    VV = zeros(nn, d);
    DD = randn(pp, d);
    for ep=1:epochs
        for t=1:nn
            z = Z(:, t);
            lambda2 = sqrt(t) * lambda2_base;
            v = OLRR_solve_v(z, DD, lambda1);
            normz = norm(z);
            u = (DD - MM)' * z / (normz * normz + 1/lambda2);
            MM = MM + z * u';
            AA = AA + v * v';
            BB = BB + z * v';
            DD = OLRR_solve_D(DD, MM, AA, BB, lambda1, lambda2);
            UU(t, :) = u';
            VV(t, :) = v';
        end
        MM = zeros(pp, d);
    end
    X = UU*VV';
    W= BuildAdjacency(thrC(X,rho));

end