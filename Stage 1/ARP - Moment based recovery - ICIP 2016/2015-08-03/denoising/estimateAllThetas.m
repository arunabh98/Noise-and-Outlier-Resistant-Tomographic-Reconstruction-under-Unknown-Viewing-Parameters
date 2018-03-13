function allThetas = estimateAllThetas(IM, allPM)
    % Finds out best estimates for all thetas given their projection
    % moments. allPM contains the projection moments as a vector
    % (including theta1, theta2)
    % 
    
    %%
    numkeep = length(allPM);
    allThetas = zeros(numkeep,1);
    possibleAngles = -180:180;
    lhs = cosd(possibleAngles)*IM(1) + sind(possibleAngles)*IM(2);
    % cosT*M10 + sinT*M01 = M(1)T. lhs is LHS of this equality. 
    % This is computed once, since it is constant for all angles
    
    
    for i = 3:numkeep
        [~, bestindex] = min(abs(lhs - allPM(i)));
        allThetas(i) = possibleAngles(bestindex);
    end
    
end
    
