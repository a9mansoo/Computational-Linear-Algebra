function [x, iter] = Jacobi(A,b, x_initial, maxiter, tol)
    % Find D-1 of A
    D = diag(diag(A));
    D_inv = inv(D);
    
    % Set the initial parameters
    iter = 0;
    x = x_initial;
    
    % 2-norm of vector r
    r = norm(b-A*x_initial);
    tolerance = tol*norm(b);
    
    % while the iter is not equal to maxiter and the residual is less 
    % than the norm(b)* tolerance
    while ( iter < maxiter && r > tolerance )
        x = x + D_inv * (b-A*x);
        r = norm(b-A*x);
        iter = iter + 1;
    end
  
end

