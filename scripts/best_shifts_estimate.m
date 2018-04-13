function estimated_projections = ...
    best_shifts_estimate(noise_amplitude, theta, image, projections)

    actual_projections = radon(image, theta);
    all_possible_shifts = -noise_amplitude:noise_amplitude;
    estimated_projections = zeros(size(projections));

    for i=1:size(projections, 2)
        min_error = inf;
        best_projection = projections(:, i);
        current_projection = actual_projections(:, i);
        for j=1:size(all_possible_shifts, 2)
            estimated_projection = ...
                circshift(projections(:, i), all_possible_shifts(j));
            error = norm(current_projection - estimated_projection);
            if error < min_error
                min_error = error;
                best_projection = estimated_projection;
            end
        end
        estimated_projections(:, i) = best_projection;
    end
end