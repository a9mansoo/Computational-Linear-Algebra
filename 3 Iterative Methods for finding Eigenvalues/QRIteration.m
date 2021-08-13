function [V, Lambda, iter] = QRIteration(A,maxiter,tol)
% initialize parameters
A0 = A;
iter = 0;
n = size(A,1);
V = eye(n);

% Loop until maxiter or if tolerance < tol
for i = 1:maxiter
   [Q,R] = qr(A0);
   
   % Eigenvectors are products of the Q's from all iterations
   V = V*Q;
   
   % Eigenvalues will be on diagonal of A0
   Lambda = diag(A0);
   
   % Update A
   A0 = R*Q;
   
   % Update iteration
   iter = iter + 1;
   
   % Check tolerance of each column of V which converge to eigenvectors and
   % respective column of Lambda which are the eigenvalues
   logical = ones(n);
    for j = 1:n
       tolerance = norm((A*V(1:n,j) - Lambda(j,1)*V(1:n,j)));
       
       if (tolerance > tol)
           logical(j) = 0;
       end
           
    end
    
    if ( all(logical) )
        break;
    end
    
end
end

