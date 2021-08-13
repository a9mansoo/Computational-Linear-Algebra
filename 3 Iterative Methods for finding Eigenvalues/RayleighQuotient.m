function [v,lambda,iter] = RayleighQuotient(A,v0,maxiter,tol)
% Initialize parameters
v = v0;
iter = 0;
lambda = transpose(v)*A*v;
tolerance = norm(A*v - lambda*v);
n = size(A,1);
disp(lambda);
while (iter < maxiter && tolerance > tol)
    
    % Solve using current lambda
    I = eye(n);
    w = (A - lambda*I) \ v;
   
    % normalize
    v = w / norm(w);
    
    %update lambda
    lambda = transpose(v)*A*v;
    

    % Update parameters and check tolerance
    iter = iter + 1;
    tolerance = norm(A*v - lambda*v);
    
end

end

