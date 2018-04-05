% Get the image.
P = imread('../images/200px-mickey.jpg');
P = imresize(P, 0.4);
P = im2double(rgb2gray(P));

% Constants.
D = dctmtx(80);
x = D*P;
x = x(:);
sigmaNoiseFraction = 0.05;
filename = ...
    '../results/moment_estimation/unknown_angles/translational_error/';
lambda  = 0.1;
rel_tol = 200;
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);

% Number of angles list.
num_theta = 30;

for o=1:1
    theta_to_write = zeros(5, num_theta(o));
    amplitude = 10;
    
    mkdir(strcat(filename, num2str(num_theta(o))));
    % Write the original image.
    imwrite(P, strcat(filename,...
        num2str(num_theta(o)), '/original_image.png'));

    % Define ground truth angles and take the tomographic projection.
    theta = datasample(0:179, num_theta(o));  
    [projections, svector] = radon(P,theta);
    theta_to_write(1, :) = theta;
 
    % Normalize s to a unit circle
    smax = max(abs(svector));
    svector = svector / smax;
    projection_length = size(projections, 1);

    % Add noise to projections.
    [projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);
    
    % Predict the angles using moment angle estimation.
    [projections, noisy_theta] = ARPord(projections, svector, sigmaNoise);
    noisy_theta = noisy_theta';

    noisy_theta = noisy_theta + theta(1) - noisy_theta(1);
    noisy_theta = process_theta(noisy_theta);
    theta_to_write(2, :) = noisy_theta;
    relative_estimated_error = norm(noisy_theta - theta)/norm(theta);

    y = projections(:);

    n = height*width;
    m = projection_length*size(noisy_theta, 2);

    % Start iteration.
    better_theta = noisy_theta;
    previous_error = inf;
    precision = 0.1;
    errors = [];

    % Reconstruct the images from projection.
    reconstructed_image = iradon(projections, noisy_theta, output_size);
    imwrite(reconstructed_image, ...
        strcat(filename, num2str(num_theta(o)), '/estimated_image.png'));

    for i=1:40
        A = radonTransform(...
            better_theta, width, height, output_size, projection_length);
        At = A';

        %run the l1-regularized least squares solver
        [reconstructed_image, status]= ...
            l1_ls(A,At,m,n,y,lambda,rel_tol,true);

        % The error we optimise.
        function_error = norm(A*reconstructed_image - y).^2 + ...
            lambda*norm(reconstructed_image, 1);
        
        % Reconstruct the image.
        reconstructed_image = reshape(reconstructed_image, [height, width]);
        reconstructed_image = D'*reconstructed_image;
        reconstructed_image(reconstructed_image < 0) = 0;

        if function_error < previous_error
            noisy_theta = better_theta;
            previous_error = function_error;

            disp(function_error);
            errors = [errors function_error];

            % Do a brute force search on all angles.
            better_theta = ...
                best_angle_alternate(precision, amplitude, noisy_theta,...
                reconstructed_image, projections);
        else
            precision = precision/1.1;
            amplitude = amplitude/1.1;

            % Do a brute force search on all angles.
            better_theta = ...
                best_angle_alternate(precision, amplitude, noisy_theta,...
                reconstructed_image, projections);
        end
    end
    
    imwrite(reconstructed_image, ...
        strcat(filename, num2str(num_theta(o)), '/reconstructed.png'));
    
    better_theta = better_theta + theta(1) - better_theta(1);
    better_theta = process_theta(better_theta);
    theta_to_write(3, :) = better_theta;
    relative_reconstructed_error = ...
        norm(better_theta - theta)/norm(theta);
    
    % Plot the function error.
    figure; plot(errors);
    saveas(gcf, ...
        strcat(filename, num2str(num_theta(o)), '/error.png'));
    
    % Write the thetas to csv file.
    theta_to_write(4, 1) = relative_estimated_error;
    theta_to_write(5, 1) = relative_reconstructed_error;
    csvwrite(strcat(filename,...
        num2str(num_theta(o)), '/thetas.csv'), theta_to_write);
end