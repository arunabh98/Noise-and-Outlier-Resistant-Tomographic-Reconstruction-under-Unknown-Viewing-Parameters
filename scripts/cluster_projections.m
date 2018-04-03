function clustered_projections = ...
    cluster_projections(projections, sigmaNoise, number_of_clusters)
    % Denoise all the projections.
    projections = denoise(projections, sigmaNoise, 117, 200);
    projections(projections < 0) = 0;
    all_projections = projections;
    
    % Cluster all the projections.
    [idx, ~] = kmeans(all_projections', number_of_clusters, 'MaxIter', 10000);    
    idx = idx'; 
    [sorted_idx, idx_order] = sort(idx);
    sorted_projections = all_projections(:, idx_order);
    
    [unique_idx, ~, ~] = unique(sorted_idx);
    number_of_projections = histc(sorted_idx, unique_idx);
    
    mean_projections = zeros(size(projections, 1), number_of_clusters);
    clustered_projections = zeros(size(projections, 1), number_of_clusters);
    c = 1;
    for i=1:size(number_of_projections, 2)
        cluster_projections = ...
            sorted_projections(:, c:c + number_of_projections(i) - 1);
        mean_projections(: , i) = mean(cluster_projections, 2);
%         clustered_projections(:, i) = denoise(mean_projections(: , i),...
%             sigmaNoise/number_of_projections(i), 10, 40);
        clustered_projections(:, i) = mean_projections(: , i);
        c = c + number_of_projections(i);
    end
    clustered_projections(clustered_projections < 0) = 0;
end