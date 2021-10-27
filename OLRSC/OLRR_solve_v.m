function v = OLRR_solve_v( z, D, lambda1)

[~, d] = size(D);
I = eye(d, d);

aux = (D' * D + 1/lambda1 * I) \ D';

v = zeros(d, 1);

eps = 1e-3;
maxIter = 1e3;

converge = false;
iter = 0;

while ~converge
    iter = iter + 1;
    
    orgv = v;
    
    v = aux * z;
    
    stopc = norm(v - orgv) / norm(v);
    
    if stopc < eps || iter > maxIter
        converge = true;
    end
end
end

