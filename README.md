# Computational-Linear-Algebra
Numerical methods associated with matrices. This repo consists of:

# 1. Method for solving linear systems:

Using LU decomposition, Cholesky Factorization and forward and backwards solves for dense and banded matrices. Also includes a method for creating a 2D Laplacian matrix. Below is an image of a meshplot showing the heat flow with 2 heat sources

![meshplot](https://user-images.githubusercontent.com/63682861/129399694-dca56b80-88a2-48b4-b1c1-7a1f53e8fc38.jpg)

# 2. Image Denoising: 

By solving linear systems derived from the discretization of a partial differential equation using iterative methods such as Gauss Siedel, Jacobi, SOR and Conjugate Gradient method.

Below is an image of denoising a 128 x 128 greyscale noisey image:

![128x128 noisey](https://user-images.githubusercontent.com/63682861/129400564-91f5bc96-118d-49d9-851b-84ecd8361355.jpg)

And the denoised image: 

![128_image](https://user-images.githubusercontent.com/63682861/129400616-420cb6a1-72df-4fcd-adc6-ae9e6ade1fb8.jpg)

# 3. Iterative Methods for solving the EigenValue Problem:

Iterative Methods known for solving the eigenvalue problem are:

- Rayleigh Quotient Iteration: finds the eigenvector and eigenvalue pair closest to the initial vector supplied.
- Power Iteration: find the eigenvector and eigenvalue pair corresponding to the largest. 
- QR Iteration: finds the set of all eigenvector and eigenvalue pairs.

A 100x100 tridiagonal matrix was supplied to the program with the initial starting column vector of [1,0,0...0]^(T) for the power iteration and [1, 1..., 1]^(T) for the Rayleigh Quotient: 

Image of eigenvector from the Rayleigh Quotient: 

![Rayleigh quotient lambda](https://user-images.githubusercontent.com/63682861/129401476-46b1a923-e9fd-4a77-bd25-fddcc659062d.jpg)

Image of eigenvector from Power Iteration: 

![Power iteration lambda](https://user-images.githubusercontent.com/63682861/129401505-92c249c1-511e-4899-b8d5-f5122feb2d2e.jpg)

Image of all eigenvectors from QR Iteration: 

![QR iteration all eigenvalues](https://user-images.githubusercontent.com/63682861/129401542-9038b065-f587-4235-a5e9-bf774dc0a629.jpg)


# 4. Specteral Clustering:

Image segmentation was performed by generating the weight matrix and forming the normalized graph Laplacian. The smallest magnitude eigenvalues and eigenvectors pairs were extracted from the normalized graph Laplacian and passed in to the kmeans clustering algorithm. 

Below is the 9 clusters found by the algorithm to identify a cell segment: 

![fig](https://user-images.githubusercontent.com/63682861/129402766-198cd2e0-7b16-4336-b9b0-e5cd4c4cb60a.jpg)







