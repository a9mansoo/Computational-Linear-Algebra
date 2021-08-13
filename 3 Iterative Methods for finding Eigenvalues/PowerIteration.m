function [v,lambda,iter] = PowerIteration(A,v0,maxiter,tol)
% Gives eigenvector with largest magnitude eigenvalue.

% Initialize parameters
iter = 0;
v = v0;

% Loop to find the eigenvector and eigenvalue
for i = 1:maxiter
    
    % Apply A, normalize and find rayleigh quotient
    w = A*v;
    v = w / norm(w);
    lambda = transpose(v)*A*v;

    % Update tolerance and iterations
    iter = iter + 1;
    tolerance = norm(A*v - lambda*v);
    
    disp(lambda);
    disp(v);
    % break out of for loop if tolerance is less than tol.
    if (tolerance < tol)
        break;
    end
    
end
    
end

