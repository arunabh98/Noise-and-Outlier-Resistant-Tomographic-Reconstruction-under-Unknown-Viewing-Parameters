modifyPMord(PMord, index, shiftEstimated)

function PMord = ...
    assemblePMord(Pgiven, kmax, shiftestimated, svector, Ord, numkeep)
    M = ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
    for k = 0:kmax
        for i = 1:size(Pgiven,2)
            shifted_proj = circshift(Pgiven(:,i), shiftestimated(i)); 
            M(i, k+1) = calculateProjectionMoment(shifted_proj, svector, k);
        end
    end
    PMord = reshape(M(:,2:(Ord+1)), Ord * numkeep, 1);
end