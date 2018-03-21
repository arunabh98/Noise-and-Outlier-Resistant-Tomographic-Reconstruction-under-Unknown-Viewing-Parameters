function denoised_projections = denoise_clusters(projections, idx) 
    for i=1:size(projections, 2)
        projections(:, i) = denoise(projections(:, i), sigmaNoise, 50, 200);
    end
end