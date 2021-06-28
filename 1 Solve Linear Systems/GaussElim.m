function x = GaussElim(A, b)

% Find A = LU
n = size(A,1);

 for k = 1:n-1
     i = k+1
         disp(i);
         mult = A(i,k) / A(k, k);
         A(i,k) = mult;
         
         for j = k+1:n
             A(i,j) = A(i,j) - mult*A(k,j);
         end
  end
 
 disp(A);
  L = tril(A,-1) + [1 0 0 0 0 ;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0; 0 0 0 0 1];
  U = triu(A);
  disp(L*U)
 % Forward solve
 z = zeros(n,1);

 for i = 1:n
     z(i) = b(i);
     for j = 1:i-1
         z(i) = z(i) - A(i,j) * z(j);
     end 
 end 
 
 % Backward solve

 x = zeros(n,1);
 
 for i = n:-1:1
     x(i) = z(i);
     for j = i+1:n
         x(i) = x(i) - A(i,j)*x(j);
     end 
     x(i) = x(i) / A(i,i);
 end
end

