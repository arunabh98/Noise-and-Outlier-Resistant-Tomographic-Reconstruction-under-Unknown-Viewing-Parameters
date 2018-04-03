% Number of angles list.
num_theta = 20000;
angle_low = 180;
angle_high = 0;

filename = strcat('../results/moment_estimation/k_means_and_unknown/',...
    num2str(num_theta), '/');

% Number of clusters.
num_clusters = [120];

for k=1:size(num_clusters, 2)
    iteration_name = ...
        strcat(num2str(num_theta), '_', num2str(num_clusters(k)));
    original_image = ...
        double(imread(strcat(filename, iteration_name, '/original_image.png')));
    estimated_image = ...
        double(imread(strcat(filename, iteration_name, '/estimated_image.png')));
    reconstructed_image = ...
        double(imread(strcat(filename, iteration_name, '/reconstructed_image.png')));
    
    error_estimated_image = inf;
    final_estimated_image = estimated_image;
    for i = angle_low:0.1:angle_high
        rotated_image = imrotate(estimated_image, i, 'crop');
        error = ...
            norm(rotated_image - original_image)/norm(original_image);
        if error < error_estimated_image
            error_estimated_image = error;
            final_estimated_image = rotated_image;
        end
    end
    disp(error_estimated_image);
    
    error_reconstructed_image = inf;
    final_reconstructed_image = reconstructed_image;
    for i = angle_low:0.1:angle_high
        rotated_image = imrotate(reconstructed_image, i, 'crop');
        error = ...
            norm(rotated_image - original_image)/norm(original_image);
        if error < error_reconstructed_image
            error_reconstructed_image = error;
            final_reconstructed_image = rotated_image;
        end
    end
    disp(error_reconstructed_image);

    figure; imshow(original_image, []);
    figure; imshow(final_estimated_image, []);
    figure; imshow(final_reconstructed_image, []);
end