function A = assembleA(thetas, order)
    %TODO: write a generalized algorithm
    
    if order == 2
        A1 = [cosd(thetas), sind(thetas), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas))];
        A2 = [zeros(size(thetas)), zeros(size(thetas)), (cosd(thetas)).^2, 2 * sind(thetas) .* cosd(thetas), (sind(thetas)).^2];
        A = [A1; A2];
    elseif order == 3
        A1 = [cosd(thetas), sind(thetas), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas))];
        A2 = [zeros(size(thetas)), zeros(size(thetas)), (cosd(thetas)).^2, 2 * sind(thetas) .* cosd(thetas), (sind(thetas)).^2, zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas))];
        A3 = [zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), (cosd(thetas)).^3, 3 * sind(thetas) .* (cosd(thetas)).^2, 3 * (sind(thetas)).^2 .* cosd(thetas), (sind(thetas)).^3];
        A = [A1; A2; A3];
    elseif order == 4
        A1 = [cosd(thetas), sind(thetas), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas))];
        A2 = [zeros(size(thetas)), zeros(size(thetas)), (cosd(thetas)).^2, 2 * sind(thetas) .* cosd(thetas), (sind(thetas)).^2, zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas))];
        A3 = [zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), (cosd(thetas)).^3, 3 * sind(thetas) .* (cosd(thetas)).^2, 3 * (sind(thetas)).^2 .* cosd(thetas), (sind(thetas)).^3, zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas))];
        A4 = [zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), zeros(size(thetas)), (cosd(thetas)).^4, 4 * sind(thetas) .* (cosd(thetas)).^3, 6 * (sind(thetas)).^2 .* (cosd(thetas)).^2, 4 * (sind(thetas)).^3 .* cosd(thetas), (sind(thetas)).^4 ];
        A = [A1; A2; A3; A4];
    end
    

    
    
end
    
