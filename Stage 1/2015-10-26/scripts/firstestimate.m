function thetas = firstestimate(actual)
    %% 1. Uniformly across entire range
    %thetas = randi([1 179], numkeep, 1);
    
    %% 2. With some uncertainty around actual values
    thetas = actual + randi([-2 2],size(actual));
    
end
