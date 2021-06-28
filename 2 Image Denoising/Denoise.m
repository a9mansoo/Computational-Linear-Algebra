% Set initial conditions

% Set outer loop
K = 8; 

% Set alpha for different image sizes
alpha_16 = 6.4 * 10^(-2);
alpha_32 = 3.2 * 10^(-2);
alpha_64 = 1.6 * 10^(-2);
%alpha_128 = 8 * 10^(-3);

% Set values for iterative methods
tol = 10^(-2);
maxiter = 20000;

% Tuned values for omega
omega_16 = 1.97;
omega_32 = 1.98;
omega_64 = 1.98;

% Create 3 random noisy images
[z_1, X_16] = set_image(16);
[z_2, X_32] = set_image(32);
[z_3, X_64] = set_image(64);
%[z_4, X_128] = set_image(128);

u_1 = FormRHS(X_16);
u_2 = FormRHS(X_32);
u_3 = FormRHS(X_64);
u_4 = FormRHS(X_128);


%[u,tot] = denoise("GS", u_4, alpha_128, omega_64, X_128);
%display_img(u);


tic
[u, tot] = denoise("Jacobi", u_1, alpha_16, omega_16, X_16);
toc
display_img(u);
disp(tot);

tic
[u, tot] = denoise("GS", u_1, alpha_16, omega_16, X_16);
toc
display_img(u);
disp(tot);

tic
[u, tot] = denoise("SOR", u_1, alpha_16, omega_16, X_16);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("CG", u_1, alpha_16, omega_16, X_16);
toc

display_img(u);
disp(tot);



tic
[u, tot] = denoise("Jacobi", u_2, alpha_32, omega_32, X_32);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("GS", u_2, alpha_32, omega_32, X_32);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("SOR", u_2, alpha_32, omega_32, X_32);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("CG", u_2, alpha_32, omega_32, X_32);
toc

display_img(u);
disp(tot);



tic
[u, tot] = denoise("Jacobi", u_3, alpha_64, omega_64, X_64);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("GS", u_3, alpha_64, omega_64, X_64);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("SOR", u_3, alpha_64, omega_64, X_64);
toc

display_img(u);
disp(tot);

tic
[u, tot] = denoise("CG", u_3, alpha_64, omega_64, X_64);
toc

display_img(u);
disp(tot);


%{
[u, tot] = denoise("Normal", u_1, alpha_16, omega_16, X_16);
display_img(u);
disp(tot);
%}

% Denoising function
function [u, tot] = denoise(iter_method, u, alpha, omega, image)
    tot = 0;
    K = 8;
    % Set values for iterative methods
    tol = 10^(-2);
    maxiter = 20000;
    if (iter_method == "Jacobi")
        for k = 0:K
            % Form matrix A(u^k)
            A = FormMatrix(u, alpha);
    
            % Form solution vector u0
            b = FormRHS(image);
    
            % Call solver to get u^k+1
            [u, iter] = Jacobi(A, b, u, maxiter, tol);
            tot = tot + iter;
        end
        
    elseif (iter_method == "GS")
        % Denoise for GS
        for k = 0:K
            % Form matrix A(u^k)
            A = FormMatrix(u, alpha);
    
            % Form solution vector u0
            b = FormRHS(image);
    
            % Call solver to get u^k+1
            [u, iter] = GS(A, b, u, maxiter, tol);
            tot = tot + iter;
        end
        
    elseif (iter_method == "SOR")
        % Denoise for SOR with tuning
        for k = 0:K
           % Form matrix A(u^k)
           A = FormMatrix(u, alpha);
    
           % Form solution vector u0
           b = FormRHS(image);
    
           % Call solver to get u^k+1
           [u, iter] = SOR(omega, A, b, u, maxiter, tol);
           tot = tot + iter;
        end
        
    elseif ( iter_method == "Normal" )
        % Denoise for \ 
        for k = 0:K
           % Form matrix A(u^k)
           A = FormMatrix(u, alpha);
    
           % Form solution vector u0
           b = FormRHS(image);
    
           % Call solver to get u^k+1
           
           u = A \ b;       
        end
        
    else
        % Denoise for CG
        for k = 0:K
          % Form matrix A(u^k)
          A = FormMatrix(u, alpha);
    
          % Form solution vector u0
          b = FormRHS(image);
    
          % Call solver to get u^k+1
          [u, iter] = CG(A, b, u, maxiter, tol);
          tot = tot + iter;
        end
    end
    
end


function [] = display_img(u)
    % Convert u back to matrix
    n = size(u,1);
    m = sqrt(n);
    img = sparse(m,m);
    for i = 1:m
        for j = 1:m
            % Calculate the row to extract
            k = i + (j-1)*m;
            img(i,j) = u(k);
        end
    end
    imagesc(img);
    colormap(gray);
end

