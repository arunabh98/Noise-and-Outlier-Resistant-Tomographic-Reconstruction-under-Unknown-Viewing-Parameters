function thetas = firstestimate(actual)
    %% 1. Uniformly across entire range
    thetas = randi([0 179], 1, size(actual, 2));
    
    %% 2. With some uncertainty around actual values
%     thetas = actual + randi([-2 2],size(actual));
    
end
