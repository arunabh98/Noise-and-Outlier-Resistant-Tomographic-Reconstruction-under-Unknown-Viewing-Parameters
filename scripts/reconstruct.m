num_theta = [10 20 30 50 100 170 500 1000];
% Get the image.
P = imread('../images/200px-mickey.jpg');
P = imresize(P, 0.4);
P = im2double(rgb2gray(P));
for o=1:8
    disp(num_theta(o));
    % Write the original image.
    imwrite(P, strcat('../results/unknown_angles/num_angles/',...
        num2str(num_theta(o)), '/original_image.png'));

    % Constants.
    snr = 1;
    num_angles = 100;
    precision = 0.1;

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
    imwrite(reconstructed_image, ...
        strcat('../results/unknown_angles/num_angles/',...
        num2str(num_theta(o)), '/estimated_image.png'));

    reconstructed_image = refine_reconstruction(projections,...
        thetasestimated', height, width, projection_length, output_size);

    imwrite(reconstructed_image, strcat('../results/unknown_angles/num_angles/',...
        num2str(num_theta(o)), '/reconstructed.png'));
end;
