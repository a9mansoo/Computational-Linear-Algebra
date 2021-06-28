function [x, iter] = CG(A,b, x_inital, maxiter, tol)

% Set initial conditions
x = x_inital;
iter = 0;
r = norm(b - A*x);
r_prev = b - A*x;
r_next = r_prev;
tolerance = tol * norm(b);

    while ( iter < maxiter && r > tolerance) 
        
        % Find beta coefficient, if 0th iteration, beta is 0 
        % Get A-orthogonal search direction pk
        if (iter == 0)
            beta = 0;
            p = r_next + beta;
        else
            beta = (transpose(r_next) * r_next)/(transpose(r_prev) * r_prev);
            p = r_next + (beta*p);
        end
        
        % get step length alpha along p
        alpha = (transpose(r_next) * r_next)/(transpose(p)*(A*p));
        
        % Update x and the residual
        x = x + alpha*p;
        r_prev = r_next;
        r_next = r_prev - (alpha*A*p);
        
        % Update looping conditions
        r = norm(r_next);
        iter = iter + 1;
    end
end

