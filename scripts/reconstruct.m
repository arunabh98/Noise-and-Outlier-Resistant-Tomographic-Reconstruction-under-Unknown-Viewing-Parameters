% Get the image.
P = imread('../images/200px-mickey.jpg');
P = imresize(P, 0.4);
% Write the original image.
imwrite(P, '../results/num_angles/30/original_image.png');
P = im2double(rgb2gray(P));
P = imresize(P, 2);

% Constants.
snr = 1;
num_angles = 100;
precision = 0.01;

% Define ground truth angles and take the tomographic projection.
possible_thetas = 0:precision:(180-precision);
theta = datasample(possible_thetas, num_angles);
projections = radon(P, theta);

% Constrain the output size of each reconstruction to be the same as the
% size of the original image, |P|.
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);
projection_length = size(projections, 1);

% Initialize angles randomly
thetasestimated = firstestimate(theta);

thetasestimated = moment_angle_estimation(projections, thetasestimated');

% Reconstruct the images from projection.
reconstructed_image = iradon(projections, thetasestimated, output_size);
reconstructed_image = imresize(reconstructed_image, 0.5);
figure; imshow(reconstructed_image);
imwrite(reconstructed_image, '../results/num_angles/30/estimated_image.png');

reconstructed_image = refine_reconstruction(projections,...
    thetasestimated, height, width, projection_length, output_size);

reconstructed_image = imresize(reconstructed_image, 0.5);
figure(); imshow(reconstructed_image);
imwrite(reconstructed_image, '../results/num_angles/30/reconstructed.png');


