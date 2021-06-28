function x = Cholesky(A, b)
n = size(A,1);
B=A;
% Find G
for k = 1:n
    A(k,k) = sqrt(A(k,k));
    
    for i = k+1:n
        A(i,k) = A(i,k)/A(k,k);
    end
    
    for j = k+1:n
        for i = j:n
            A(i,j) = A(i,j) - A(i,k) * A(j,k);
        end
    end    
end

% Get G and GT

G = tril(A);
G_trans = transpose(G);

z=zeros(n,1);

% Forward solve
 for i = 1:n
     z(i) = b(i);
     for j = 1:i-1
         z(i) = z(i) - G(i,j) * z(j);
     end 
     z(i) = z(i) / G(i,i);
 end
 
 x=zeros(n,1);
 
% Backward solve
for i = n:-1:1
     x(i) = z(i);
     for j = i+1:n
         x(i) = x(i) - G_trans(i,j)*x(j);
     end 
     x(i) = x(i) / G_trans(i,i);
end
end

