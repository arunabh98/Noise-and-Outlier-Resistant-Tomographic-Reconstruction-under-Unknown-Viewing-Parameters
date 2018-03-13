
candidate_index = besttheta2candidate_index;

expectedanswer = expectedanswer_original;
expectedanswer(2) = expectedanswer_original(candidate_index);
expectedanswer(candidate_index) = expectedanswer_original(2);

M = M_original;
M(2, :) = M_original(candidate_index, :);
M(candidate_index, :) = M_original(2, :);

%function allthetas = EstimateTheta2(M)

theta1 = 0;
bestanswer = zeros(numkeep);
besterror = Inf;

bestIMreconstructed = zeros(2,1);
bestIMassumed = zeros(2,1); % IMassumed in the case where *reconstruction* is best. 
IMassumed = zeros(2,1);

rangethetas = 1:179;
ctr = 1;
totalerrorvector = zeros(size(rangethetas));
%Rearranged thetas to put appropriate angle in position (2)

PM1 = M(:, 2); % 1 indicates order 

for theta2 = rangethetas    
    %%  Calculate (first order) image moments from 2 projection moments and angles

    A = [cosd(theta1), sind(theta1); cosd(theta2), sind(theta2)];
    IMassumed = A \ PM1(1:2);   % = A^-1 * PM, using first two values

    %%  Estimate all thetas from image moments and remaining projection moments

    estimatedthetas =  estimateAllThetas(IMassumed, PM1);
    estimatedthetas(1) = 0;
    estimatedthetas(2) = theta2; 

    [IMreconstructed, error] = evaluateGoodnessOfThetas(estimatedthetas, IMassumed, PM1);

    totalerrorvector(ctr) = error;
    ctr = ctr + 1;
    %% Evaluate error in moment estimates. Return thetas with lowest error
    if error < besterror
        bestanswer = estimatedthetas;
        besterror = error;
        bestIMreconstructed = IMreconstructed;
        bestIMassumed = IMassumed;
    end

end

clc;
disp('Expected | Computed');
[expectedanswer ,estimatedthetas]
besterror