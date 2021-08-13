%
% CS475/675: Assignment 4
%
%   Cell image segmentation
%


%
% Read in a block from a cell image
%
U = imread('cellimage.tif');
U = U(90:190,190:290);
U = double(U);
U = U/max(U(:));


%
% Create the normalized graph Laplacian from image U
% You will need to implement this function for ***part (a).***

NL = CreateImageGraph(U);

disp(NL);
%
% Perform normalized spectral clustering
%

%*** Provide your code here for part (b)! ***
K = 9;

[V, D] = eigs(NL, K, 'smallestabs');

n = size(V,1);

for i = 1:n
   V(i,:) = V(i,:)/norm(V(i,:)); 
end
    


index = kmeans(V, K, 'Replicates', 20);

%result should be the variable 'index' produced by Matlabs kmeans command, 
%i.e., a vector of length m*n containing the cluster index for each pixel



%
% Extract segments for the expected cell region in a simple way
%
Clusters = reshape(index,size(U,1),size(U,2));

Disk = fspecial('disk',floor(size(U,1)/2));
Disk = Disk>0;

Cell = zeros(size(U));
for k=1:K
    seg_size = nnz(Clusters==k);
    overlap = (Clusters==(Disk*k));
    in_size = nnz(overlap);
    if in_size == seg_size,
        Cell = Cell + (Clusters==k);
    end
end
Cell = 2*(Cell-0.5);


%
% Visualize segmentation results
%
figure(1);

%input image
subplot(1,3,1);
imshow(U,[]);

%generated clusters
subplot(1,3,2);
imshow(Clusters,[]);

%segmented result
subplot(1,3,3);
imshow(U,[]);
hold on;
contour(Cell,[0 0],'r', 'linewidth', 1.5);
hold off;
