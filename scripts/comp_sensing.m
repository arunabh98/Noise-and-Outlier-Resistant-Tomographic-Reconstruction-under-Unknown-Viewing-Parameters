% Get the image.
P = phantom(128);
D = dctmtx(128);
x = D*P;
x = x(:);
snr = 1;

imwrite(P, 'results/angles_starved/18/original_image.png');

% Define ground truth angles and take the tomographic projection.
theta = 0:10:179;  
projections = radon(P,theta);
% Add noise to projections. 
projections = awgn(projections, snr);
y = projections(:);

% Constrain the output size of each reconstruction to be the same as the
% size of the original image, |P|.
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);
projection_length = size(projections, 1);

% Add noise to the thetas of some fixed amplitude.
amplitude = 2;

% Uniform continious noise.
noisy_theta = theta + amplitude*rand(1, size(theta, 2)) - amplitude/2;

n = height*width;
m = projection_length*size(noisy_theta, 2);

lambda  = 0.1; % regularization parameter
rel_tol = 200; % relative target duality gap

% Start iteration.
better_theta = noisy_theta;
previous_error = inf;
precision = 0.1;
errors = [];

for i=1:40
    A = radonTransform(...
        better_theta, width, height, output_size, projection_length);
    At = A';
    
    tic;
    %run the l1-regularized least squares solver
    [reconstructed_image, status]=l1_ls(A,At,m,n,y,lambda,rel_tol,true);
    toc;
    
    % The error we optimise.
    function_error = norm(A*reconstructed_image - y).^2 + ...
        lambda*norm(reconstructed_image, 1);
    
    reconstructed_image = reshape(reconstructed_image, [height, width]);
    reconstructed_image = D'*reconstructed_image;
    
    current_error = norm(reconstructed_image - P);
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

% Show the reconstructed image.
figure; imshow(reconstructed_image);

figure; plot(errors);
saveas(gcf,'results/angles_starved/18/error.png');

imwrite(reconstructed_image, 'results/angles_starved/18/reconstructed_image.png');



