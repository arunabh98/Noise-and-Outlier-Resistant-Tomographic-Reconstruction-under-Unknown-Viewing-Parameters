% Get the image.
P = imread('../images/200px-mickey.jpg');
P = imresize(P, 0.4);
P = im2double(rgb2gray(P));

% Constants.
D = dctmtx(80);
x = D*P;
x = x(:);
sigmaNoiseFraction = 0.5;
filename = ...
    '../results/moment_estimation/k_means_and_unknown/';
lambda  = 0.1;
rel_tol = 200;
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);
iteration_name = '500_60/';

% Number of angles list.
num_theta = 3000;

for o=1:1
    amplitude = 10;

    % Define ground truth angles and take the tomographic projection.
%     theta = datasample(0:179, num_theta(o));  
    theta = 1:0.36:179;
    [projections, svector] = radon(P,theta);
 
    % Normalize s to a unit circle
    smax = max(abs(svector));
    svector = svector / smax;
    projection_length = size(projections, 1);

    % Add noise to projections.
    [all_projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);
    
    % Cluster the projections into groups.
    % Using k-means clustering and denoise them.
    projections = cluster_projections(all_projections, sigmaNoise);
    
    % Predict the angles using moment angle estimation.
    noisy_theta = ARPord_Kmeans(projections, svector);
    noisy_theta = noisy_theta';

    noisy_theta = noisy_theta + theta(1) - noisy_theta(1);
    noisy_theta = process_theta(noisy_theta);

    y = projections(:);

    n = height*width;
    m = projection_length*size(noisy_theta, 2);

    % Start iteration.
    better_theta = noisy_theta;
    previous_error = inf;
    precision = 0.1;
    errors = [];

    % Reconstruct the images from projection.
    estimated_image = iradon(projections, noisy_theta, output_size);

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
    
    better_theta = better_theta + theta(1) - better_theta(1);
    better_theta = process_theta(better_theta);
    
    % Store all the results.
    imwrite(reconstructed_image, strcat(filename,...
        iteration_name, '/reconstructed_image.png'));
    imwrite(estimated_image, strcat(filename,...
        iteration_name, '/estimated_image.png'));
    imwrite(P, strcat(filename,...
        iteration_name, '/original_image.png'));
    fileID = fopen(strcat(filename, iteration_name,'thetas.txt'), 'w');
    fprintf(fileID,'%16s %20s\n','estimated angles', 'reconstructed angles');
    noisy_theta = noisy_theta';
    better_theta = better_theta';
    fprintf(fileID,'%3.3s %3.3s\n',[noisy_theta better_theta]);
    fclose(fileID);
end;