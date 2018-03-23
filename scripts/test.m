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
    '../results/moment_estimation/k_means_and_unknown/';
lambda  = 0.1;
rel_tol = 200;
output_size = max(size(P));
height = size(P, 1);
width = size(P, 2);

% Number of angles list.
num_theta = 100;
    amplitude = 10;

    % Define ground truth angles and take the tomographic projection.
    theta = 1:180;  
    [projections, svector] = radon(P,theta);
 
    % Normalize s to a unit circle
    smax = max(abs(svector));
    svector = svector / smax;
    projection_length = size(projections, 1);

    % Add noise to projections.
    [all_projections, sigmaNoise] = add_noise(projections, sigmaNoiseFraction);
    
    % Cluster the projections into groups.
    % Using k-means clustering and denoise them.
    number_of_clusters = 30;
    % Cluster all the projections.
    [idx, ~] = kmeans(all_projections', number_of_clusters);    
    idx = idx'; 
%     projections = cluster_projections(all_projections, sigmaNoise);
%     
%     % Predict the angles using moment angle estimation.
%     noisy_theta = ARPord_Kmeans(projections, svector);
%     noisy_theta = noisy_theta';
% 
%     noisy_theta = noisy_theta + theta(1) - noisy_theta(1);
%     noisy_theta = process_theta(noisy_theta);