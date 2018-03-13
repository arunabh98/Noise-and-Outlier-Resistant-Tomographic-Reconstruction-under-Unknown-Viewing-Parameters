
function A = assembleA(thetas, k)
    % thetas in degrees), k^th image moment to be calculated
    thetas = thetas(1:k+1); % The rest can be discarded
    
    A = zeros(k + 1, k + 1);
    for l = 0:k
        A(:,l+1) = nchoosek(k, k-l) * (cosd(thetas)).^(k-l) .* (sind(thetas)).^l;
        % l+1 because 1-indexed matrix
    end
    
end