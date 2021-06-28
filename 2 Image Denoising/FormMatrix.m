function [A] = FormMatrix(u,alpha)

% Initial parameters
n = size(u,1); 

m = sqrt(n);

h = 1 /(m+1);

A = sparse(n,n);

beta = 10^(-6);

for i = 1:m
   for j = 1:m
       
       % Calculate row for A
       row = i + (j-1)*m;
       % Calculate column for coefficients
       
       AC_col = i + (j-1)*m;
       AW_col = (i-1) + (j-1)*m;
       AE_col = (i+1) + (j-1)*m;
       AS_col = i + (j-1-1)*m;
       AN_col = i + (j+1-1)*m;
       
       % Calculate AW
       % u i,j
       u_1row = i + (j-1)*m;
       
       u_1 = 0;
       u_2 = 0;
       u_3 = 0;
       u_4 = 0;
       u_5 = 0;
       u_6 = 0;
       u_7 = 0;       
       
       u_1 = u(u_1row);
       % check if these points are on the boundary
       if ( (i-1) > 0 )
           % u i-1, j
           u_2row = (i-1) + (j-1)*m;
           u_2 = u(u_2row);
       end
       
       if ( (j-1) > 0 )
           % u i, j-1
           u_3row = i + (j-1-1)*m;
           u_3 = u(u_3row);
       end
       
       if ( (i-1) > 0 && (j+1) < (m+1) )
           % u i-1,j+1
           u_4row = (i-1) + (j+1-1)*m;
           u_4 = u(u_4row);
       end
       
       if ( (i+1) < (m+1) )
           % u i+1, j
           u_5row = (i+1) + (j-1)*m;
           u_5 = u(u_5row);
       end
       
       if ( ((i+1) < (m+1)) && ((j-1) > 0) )
           % u i+1,j-1
           u_6row = (i+1) + (j-1-1)*m;
           u_6 = u(u_6row);
       end
       
       if ( (j+1) < (m+1) )
           % u i, j+1
           u_7row = i + (j+1-1)*m;
           u_7 = u(u_7row);
       end
       
       % Calculate AW
       sum_1 = (u_1 - u_2);
       sum_2 = (u_1 - u_3);
       sum_3 = (u_1 - u_2);
       sum_4 = (u_4 - u_2);
       sum_5 = (u_5 - u_1);
       sum_6 = (u_5 - u_6);
       sum_7 = (u_7 - u_1);
       sum_8 = (u_6 - u_3);
       sum_9 = (u_7 - u_4);
       const = -1 * (alpha / (h^2));
       
       lhs = 1 / (sqrt( (sum_1/h)^2 + (sum_2/h)^2 + beta ));
       rhs = 1 / (sqrt( (sum_3/h)^2 + (sum_4/h)^2 + beta ));
       AW = const * (lhs + rhs);
      
       % Calculate AE
       lhs = 1 / (sqrt( (sum_5/h)^2 + (sum_6/h)^2 + beta ));
       rhs = 1 / (sqrt( (sum_5/h)^2 + (sum_7/h)^2 + beta ));
       AE = const * (lhs + rhs);
       
       % Calculate AS
       lhs = 1 / (sqrt( (sum_1/h)^2 + (sum_2/h)^2 + beta ));
       rhs = 1 / (sqrt( (sum_8/h)^2 + (sum_2/h)^2 + beta ));
       AS = const * (lhs + rhs);
       
       % Calculate AN
       lhs = 1 / (sqrt( (sum_5/h)^2 + (sum_7/h)^2 + beta ));
       rhs = 1 / (sqrt( (sum_9/h)^2 + (sum_7/h)^2 + beta ));
       AN = const * (lhs + rhs);
       
       % Calculate AC 
       AC = (-1 * (AW + AE + AS + AN)) + 1;
       
       % Check if indices are on boundary
       if ( (i-1) > 0 )
           A(row, AW_col) = AW;
       end
       
       if ( (i+1) < (m+1) )
           A(row, AE_col) = AE;
       end
       
       if ( (j-1) > 0 )
           A(row, AS_col) = AS;
       end
       
       if ( (j+1) < (m+1) )
           A(row, AN_col) = AN;
       end
       
       A(row, AC_col) = AC; 
       
   end
  
end

end

