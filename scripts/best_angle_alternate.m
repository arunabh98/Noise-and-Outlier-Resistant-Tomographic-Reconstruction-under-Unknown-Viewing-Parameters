function best_angle = ...
    best_angle_alternate(precision, noise_amplitude, noisy_theta,...
    image, projections)
    
    % Generate all possible angles.
    all_possible_angles = ...
        zeros(size(noisy_theta, 2), ceil(noise_amplitude/precision) + 1);
    parfor i=1:size(noisy_theta, 2)
        possible_angles = ...
            zeros(1, ceil(noise_amplitude/precision) + 1);
        dummy = ...
            noisy_theta(i)-noise_amplitude/2:precision:noisy_theta(i)+noise_amplitude/2;
        possible_angles(1:size(dummy, 2)) = ...
            noisy_theta(i)-noise_amplitude/2:precision:noisy_theta(i)+noise_amplitude/2;
        all_possible_angles(i, :) = possible_angles;
    end
    
    % Start choosing the best angle.
    best_angle = all_possible_angles(:, 2)';
    
    parfor i=1:size(all_possible_angles, 1)
        min_error = inf;
        pos_best_theta = 1;
        for j=1:size(all_possible_angles, 2)
            best_angle(i) = all_possible_angles(i, j);
            estimated_projections = radon(image, all_possible_angles(i, j));
            error = norm(projections(:, i) - estimated_projections);
            if error < min_error
                min_error = error;
                pos_best_theta = j;
            end
        end
        best_angle(i) = all_possible_angles(i, pos_best_theta);
    end
end
