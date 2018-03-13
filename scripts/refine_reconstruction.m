function [reconstructed_image, better_theta] = refine_reconstruction(projections,...
    noisy_theta, height, width, projection_length, output_size)

    y = projections(:);
    D = dctmtx(height);
    amplitude = 12;

    lambda  = 0.1; % regularization parameter
    rel_tol = 200; % relative target duality gap
    
    n = height*width;
    m = projection_length*size(noisy_theta, 2);

    % Start iteration.
    better_theta = noisy_theta;
    previous_error = inf;
    precision = 0.1;

    for i=1:40
        A = radonTransform(...
            better_theta, width, height, output_size, projection_length);
        At = A';

        %run the l1-regularized least squares solver
        [reconstructed_image, ~]=l1_ls(A,At,m,n,y,lambda,rel_tol,true);

        % The error we optimise.
        function_error = norm(A*reconstructed_image - y).^2 + ...
            lambda*norm(reconstructed_image, 1);

        reconstructed_image = reshape(reconstructed_image, [height, width]);
        reconstructed_image = D'*reconstructed_image;

        if function_error < previous_error
            noisy_theta = better_theta;
            previous_error = function_error;

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
end