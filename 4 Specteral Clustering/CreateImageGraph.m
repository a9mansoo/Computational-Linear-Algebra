function [NL] = CreateImageGraph(U)
    % Set parameter values
    m = size(U,1);
    n = size(U,2);
    
    sigma_dist = 150;
    sigma_intensity = 0.002;
    
    % size of the W and D matrix
    length = m*n;
    u = sparse(length,1);
    
    % Flatten U
    for i = 1:m
        for j = 1:n
            % row-wise flatten
            k = i + (j-1)*m;
            u(k,1) = U(i,j);
        end    
    end
    
    W = sparse(length, length);
    
    % Calculate weight matrix W 
    for i = 1:m
        for j = 1:n
            
            % Calculate row
            k = i + (j-1)*m;
           
             
          % Check neighbourhood of current pixel
          % u i+1 j
         
          if ( (i+1) < (m+1) )
             
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i+1, j] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = (i+1) + (j-1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
              
          end
          
          % u i-1 j
       
          if ( (i-1) > 0 )    
              
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i-1, j] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = (i-1) + (j-1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
              
          end
          
         
          % u i-1 j-1
          
          if ( (i-1) > 0 && (j-1) > 0 )
              
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i-1, j-1] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = (i-1) + (j-1-1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
             
          end
          
          % u i  j-1
          
          if ( (j-1) > 0 )
             
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i, j-1] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = i + (j-1-1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
             
          end
          
          % u i+1 j-1
          
          if ( (i+1) < (m+1) && (j-1) > 0 )
             
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i+1, j-1] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = (i+1) + (j-1-1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
             
          end
          
          % u i-1 j+1
          
          if ( (i-1) > 0 && (j+1) < (n+1) )
              
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i-1, j+1] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = (i-1) + (j-1+1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
           
          end
         
          % u i j+1
          
          if ( (j+1) < (n+1) )
              
             % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i, j+1] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = i + (j-1+1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2)); 
              
          end
         
          % u i+1 j+1
          
          if ( (i+1) < (m+1) && (j+1) < (n+1) )
              
              % Calculate weight 
              exp_1 = -1 * ((norm( [i,j] - [i+1, j+1] ))^(2)) / sigma_dist;
              row_1 = i + (j-1)*m;
              row_2 = (i+1) + (j-1+1)*m;
              exp_2 = -1 * ((norm ( u(row_1,1) - u(row_2,1) ) )^(2)) / sigma_intensity;
              W(k, row_2) = (exp(exp_1)) * (exp(exp_2));
          
          end
          
        end
    end
    
    
    % Calculate degree D based on weight matrix W
    D = sparse(length, length);
    
    for i = 1:length
        D(i,i) = sum(W(i,:));
    end
   
    
    % Calculate normalized Laplacian 
    NL = speye(length) - (D^(-1/2) * W * (D^(-1/2)));

end

