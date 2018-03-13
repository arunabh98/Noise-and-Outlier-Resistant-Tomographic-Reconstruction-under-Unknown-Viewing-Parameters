function thetas = firstestimate(actual, deviation)
    %% 1. Uniformly across entire range
    thetas = randi([0 179], 1, size(actual, 2));
%     reverse_thetas = fliplr(thetas);
%     reverse_thetas = 180 - reverse_thetas;
%     thetas = [thetas reverse_thetas];
    %% 2. With some uncertainty around actual values
    thetas = actual + randi([-deviation deviation],size(actual));
    
end
