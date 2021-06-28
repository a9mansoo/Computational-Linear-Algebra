function [A, b] = Lap2D(m)

n = m^2;

A = zeros(n,n);

% Ti,j Ti-1,j Ti+1,j Ti,j-1, Ti,j+1
temp = [4 -1 -1 -1 -1];

for i = 1:m
    for j = 1:m        
        % Find row
        k = (j-1)*m + i;
        
        % Find the columns to place temp into
        T1 = (j-1)*m + i;
        T2 = (j-1)*m + (i-1);
        T3 = (j-1)*m + (i+1);
        T4 = (j-1-1)*m + i;
        T5 = (j+1-1)*m + i;
        
        A(k, T1) = temp(1,1);
        
        if (i-1 > 0)
            A(k, T2) = temp(1,2);
        end
        
        if (i + 1 < m+1)
            A(k, T3) = temp(1,3);
        end
        
        if ( j-1 > 0)
            A(k, T4) = temp(1,4);
        end
        
        if ( j+ 1 < m + 1)
            A(k, T5) = temp(1,5);
        end
    end
end

h = 1/(m+1);

A = 1/h^(2) * A;

% Calculate b vector

b = zeros(n,1);
point_1 = [0.35 0.6];
point_2 = [0.8 0.25];

for i = 0:m+1
    for j = 0:m+1
        
        % Row in f_i,j
        k = (j-1)*m + i;
        
        % Calculate point
        x_i = i*h;
        y_i = j*h;
        
        eval_1 = [x_i y_i] - point_1;
        eval_2 = [x_i y_i] - point_2;
        
        % Calculate disatance
        if (norm(eval_1) <= 0.1 || norm(eval_2) <= 0.1)
            b(k) = 1;
        end
    
    end
end

end
