num_theta = [20 30 50 80 100];
% Get the image.
P = imread('../images/200px-mickey.jpg');
P = imresize(P, 0.4);
P = im2double(rgb2gray(P));
parfor o=3:3
    theta_to_write = zeros(3, num_theta(o));
    disp(num_theta(o));
    % Write the original image.
    imwrite(P, strcat('../results/moment_estimation/unknown_angles/5_percent_noise/',...
        num2str(num_theta(o)), '/original_image.png'));

    % Constants.
    snr = 1;
    num_angles = num_theta(o);
    precision = 1;
    sigmaNoiseFraction = 0.05;

    % Define ground truth angles and take the tomographic projection.
    possible_thetas = 0:precision:(180-precision);
    theta = datasample(possible_thetas, num_angles);
    [projections, svector]  = radon(P, theta);
    theta_to_write(1, :) = theta;
    
    % Add noise
    ref = std(projections(:));
    sigmaNoise = sigmaNoiseFraction * ref;
    noise = normrnd(0, sigmaNoise, size(projections));
    projections = projections + noise;
    
    % Normalize s to a unit circle
    smax = max(abs(svector));
    svector = svector / smax;

    % Constrain the output size of each reconstruction to be the same as the
    % size of the original image, |P|.
    output_size = max(size(P));
    height = size(P, 1);
    width = size(P, 2);
    projection_length = size(projections, 1);

    [projections, thetasestimated] = ARPord(projections, svector, sigmaNoise);
    theta_to_write(2, :) = thetasestimated;

    % Reconstruct the images from projection.
    reconstructed_image = iradon(projections, thetasestimated, output_size);
    imwrite(reconstructed_image, ...
        strcat('../results/moment_estimation/unknown_angles/5_percent_noise/',...
        num2str(num_theta(o)), '/estimated_image.png'));

    [reconstructed_image, better_theta] = refine_reconstruction(projections,...
        thetasestimated', height, width, projection_length, output_size);
    theta_to_write(3, :) = better_theta;
    imwrite(reconstructed_image, strcat('../results/moment_estimation/unknown_angles/5_percent_noise/',...
        num2str(num_theta(o)), '/reconstructed.png'));
   
    csvwrite(strcat('../results/moment_estimation/unknown_angles/5_percent_noise/',...
        num2str(num_theta(o)), '/thetas.csv'), theta_to_write);
end;
