function [u0] = FormRHS(X)

m = size(X,1);

n = m^2;

u0 = sparse(n,1);

for i = 1:m
    for j = 1:m
    
        % Find row to place it in u0
        k = i + (j-1)*m;
        u0(k) = X(i,j);
    end 
end




end

