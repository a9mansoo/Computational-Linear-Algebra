% Create a tridiagonal matrix for n = 100
n = 10;
A = zeros(n,n);

for i = 1:n
    A(i,i) = 2;
end 

% diag(v,k) places the elements of vector v on the kth diagonal
v = -1*ones(1,n-1);
lower = diag(v,-1);
lower_2 = diag(v, -2);
upper = diag(v,1);
upper_2 = diag(v,2);
lower_3 = diag(v,-3);
upper_3 = diag(v,3);

A = A + lower + upper + upper_2(1:n,1:n) + lower_2(1:n,1:n) + lower_3(1:n,1:n) + upper_3(1:n,1:n);

e1 = zeros(9,1);
e1(1,1) = 1;

disp(A);
v = A(2:n,2);
v = v + sign(v(1,1)) * norm(v)*e1;
Q = eye(9) - (2*v*transpose(v)) / (transpose(v)*v);

disp(Q);

[E,V] = eigs(A);

J = eye(n) - 0.5*A; 

[E_2, V_2] = eigs(J); 

disp(E);
disp(V);

disp(E_2);
disp(V_2);

% Set maxiter and tolerance
maxiter = 10000;
tol = 10^(-4);



% i) Power iteration
v0 = zeros(n,1);
v0(1,1) = 1;
[v_1, lambda_1, iter_1] = PowerIteration(A, v0, maxiter, tol);


% Plot for PowerIteration 
plot(v_1);
title(['Power Iteration ',' lambda ', num2str(lambda_1),' ', ' iterations ', num2str(iter_1)]);



% ii) Rayleigh Quotient
v0 = ones(n,1);
[v_2, lambda_2, iter_2] = RayleighQuotient(A, v0, maxiter, tol);

% Plot for RayleighQuotient
plot(v_2);
title(['Rayleigh Quotient ',' lambda ', num2str(lambda_2),' ', ' iterations ', num2str(iter_2)]);

% iii) QR Iteration 
[V, Lambda, iter_3] = QRIteration(A, maxiter, tol);


% Plot for all columns of V from QR Iteration
asc_lam = flipud(Lambda);
asc_v = fliplr(V);

plot(asc_v);
title('Eigenvalues of A');

% Plot v(20), v(40), v(60), v(80)
v_20 = asc_v(1:n,20);

v_40 = asc_v(1:n,40);

v_60 = asc_v(1:n,60);

v_80 = asc_v(1:n,80);


% Plot for v(20)
plot(v_20);
title(['eigenvector 20 ','lambda ', num2str(asc_lam(20,1)), ' Iteration ', num2str(iter_3)])

% Plot for v(40)
plot(v_40);
title(['eigenvector 40 ', 'lambda ', num2str(asc_lam(40,1)), ' Iteration ', num2str(iter_3)])

% Plot for v(60)
plot(v_60);
title(['eigenvector 60 ', 'lambda ', num2str(asc_lam(60,1)), ' Iteration ', num2str(iter_3)])

% Plot for v(80)
plot(v_80);
title(['eigenvector 80 ', 'lambda ', num2str(asc_lam(80,1)), ' Iteration ', num2str(iter_3)])


