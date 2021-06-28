function x = BandGE(A, b, p, q)
    % Banded LU factorization
    n = size(A,1);

    for k = 1:n-1
        for i = k+1:min(k+p,n)
            A(i,k) = A(i,k)/A(k,k);
        end
        for i = k+1:min(k+p,n)
            for j = k+1:min(k+q,n)
                A(i,j) = A(i,j) - A(i,k)*A(k,j);
            end
        end
    end
    
      
    % Forward solve
    z=zeros(n,1);
    for i = 1:n
        z(i) = b(i);
        for j = max(1, i-p):i-1
            z(i) = z(i) - A(i,j) * z(j);
        end 
    end
    
    
    x = zeros(n,1);
    % Backward solve
    for i = n:-1:1
     x(i) = z(i);
     for j = i+1:min(i+q,n)
         x(i) = x(i) - A(i,j)*x(j);
     end 
     x(i) = x(i) / A(i,i);
    end
end

