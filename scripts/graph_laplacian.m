% Get the image.
% P = phantom(128);
P = imread('200px-mickey.jpg');
P = im2double(rgb2gray(P));
P = imresize(P, 2);

% Constants
output_size = max(size(P));
d = 2;
epsilon = 6000;
precision = 0.01;
num_angles = 500;

% Define ground truth angles and take the tomographic projection.
possible_thetas = 0:precision:(180-precision);
theta = datasample(possible_thetas, num_angles);
projections = radon(P, theta);
projections = [projections flipud(projections)];

% Compute the norm 2 of each projection with each other. 
projections2 = sum(projections.^2,1);
distance = ...
    repmat(projections2, size(projections, 2), 1) + ...
    repmat(projections2', 1, size(projections, 2)) - ...
    2*(projections')*projections;

% Compute Graph Laplacian Matrix.
W = exp(-distance/(2*epsilon));
D = diag(sum(W, 2));
W = D\W/D;
D = diag(sum(W, 2));
L = D\W - eye(size(theta, 2)*2);

% Compute the eigenvalue and reduce the dimension.
options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(L,d+1,0,options);
Y = Y(:,2:d+1)'*sqrt(size(projections, 2));

% Compute the angles in the 2 dimensional space.
A = sum(Y < 0);
B = sum((A == 2), 1);
B(B == 1) = 2;
B = -1.*(B - 1);
angles = B.*(sign(Y(1, :)).*0.5.*pi.*(sign(-Y(2, :)) + 1) + atan(Y(1, :)./Y(2, :)));

% Sort the projections.
[a_sorted, a_order] = sort(angles);
new_projections = projections(:, a_order);

% Reconstruct the image.
pred_angles = 0:(360/(2*num_angles)):360;
pred_angles = pred_angles(1:size(angles, 2));
reconstructed_image = iradon(new_projections, pred_angles, output_size);
reconstructed_image = imresize(reconstructed_image, 0.5);
imshow(reconstructed_image);


