[A_8, b_8] = Lap2D(8);
[A_16, b_16] = Lap2D(16);
[A_24, b_24] = Lap2D(24);
[A_32, b_32] = Lap2D(32);

[lower_8, upper_8] = bandwidth(A_8);
[lower_16, upper_16] = bandwidth(A_16);
[lower_24, upper_24] = bandwidth(A_24);
[lower_32, upper_32] = bandwidth(A_32);

% Gaussian Elimination
tic
x_8_gauss = GaussElim(A_8, b_8);
toc

tic
x_16_gauss = GaussElim(A_16, b_16);
toc

tic
x_24_gauss = GaussElim(A_24, b_24);
toc

tic
x_32_gauss = GaussElim(A_32, b_32);
toc

% Cholesky factorization
tic
x_8_chol = Cholesky(A_8, b_8);
toc

tic
x_16_chol = Cholesky(A_16, b_16);
toc

tic
x_24_chol = Cholesky(A_24, b_24);
toc

tic
x_32_chol = Cholesky(A_32, b_32);
toc

% BandGE
tic
x_8_band = BandGE(A_8, b_8, lower_8, upper_8);
toc

tic
x_16_band = BandGE(A_16, b_16, lower_16, upper_16);
toc

tic
x_24_band = BandGE(A_24, b_24, lower_24, upper_24);
toc

tic
x_32_band = BandGE(A_32, b_32, lower_32, upper_32);
toc