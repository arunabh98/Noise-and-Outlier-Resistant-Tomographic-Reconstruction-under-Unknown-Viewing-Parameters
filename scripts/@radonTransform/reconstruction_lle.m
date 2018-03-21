% Get the image.
P = phantom(256);
P = imresize(P, 5);

% Constants
output_size = max(size(P));
K = 10;
d = 2;
num_angles = 500;
precision = 0.01;

% Define ground truth angles and take the tomographic projection.
possible_thetas = 0:precision:(180-precision);
theta = datasample(possible_thetas, num_angles);
projections = radon(P, theta);
projections = [projections flipud(projections)];

projections2 = sum(projections.^2,1);
distance = ...
    repmat(projections2, size(projections, 2), 1) + ...
    repmat(projections2', 1, size(projections, 2)) - ...
    2*(projections')*projections;
[sorted,index] = sort(distance);
neighborhood = index(2:(1+K),:);
neighbour_indices = neighborhood';

W = zeros(size(projections, 2), size(projections, 2));
for i=1:size(projections, 2)
    Z = ...
        projections - repmat(projections(:, i), [1, size(projections, 2)]);
    Z(:, i) = [];
    C = Z'*Z;
    [w, ~] = linsolve(C, ones(size(projections, 2) - 1, 1));
    W(i, i) = 0;
    W(1:i-1, i) = w(1:i-1);
    W(i+1:end, i) = w(i:end);
    multi = zeros(size(projections, 2), 1);
    multi(neighbour_indices(i, :), 1) = 1;
    W(:, i) = W(:, i).*multi;
    W(:, i) = W(:, i)/sum(W(:, i));
end

W = sparse(W);
M = sparse((eye(size(projections, 2)) - W)'*(eye(size(projections, 2)) - W));
options.disp = 0; options.isreal = 1; options.issym = 1; 
[Y,eigenvals] = eigs(M,d+1,0,options);
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
pred_angles = pred_angles(1, size(angles, 2));
reconstructed_image = iradon(new_projections, pred_angles, output_size);
reconstructed_image = imresize(reconstructed_image, 0.2);
imshow(reconstructed_image);





