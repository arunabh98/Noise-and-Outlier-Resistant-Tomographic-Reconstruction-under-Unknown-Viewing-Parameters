%%
% Reconstruction from projection data from unknown angles
% In this script we'll implement angle recovery. From there, reconstruction
% is quite directly possible

clear; clc;
rng(100);

%% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
thetas = (0:179); % Projections to sample from
numkeep = 31; % Number of angles to recompute from
sigmaNoiseFraction = 0.001;


%% %%%%%%%%%%%%% A. GENERATE INPUT DATA %%%%%%%%%%%%%%%%%

%%%% Define image I %%%%
%I = phantom(200);
%I = rgb2gray(imread('../images/256px-Earthrise.jpg'));
%I = rgb2gray(imread('../images/200px-beach.jpg'));
I = rgb2gray(imread('../images/200px-woody.jpg'));
% figure();
% imshow(I,'InitialMagnification','fit');



%%%% Calculate projections %%%%

[P, svector] = radon(I, thetas); % P(s, theta)



% figure();
% imshow(P / max(max(P)), [], 'Xdata', [0 179] , 'Ydata', svector, 'InitialMagnification','fit');
% iptsetpref('ImshowAxesVisible','on')
% xlabel('Angle \theta')
% ylabel('s')
% title('Projections P_\theta(s)')

% Normalize s to a unit circle
smax = max(abs(svector));
svector = svector / smax;

% sample uniformly
%keep = randperm(size(P,2), numkeep);

%sample from a skewed distribution
weights = [0.1 0.5 0.07 0.3 0.03]; 
partitionsize = length(thetas)/length(weights);

keep = [];
for i = 1:length(weights)
    start = partitionsize * (i - 1);
    numsample = round(numkeep * weights(i));
    keep_partition = start + randperm(partitionsize, numsample);
    keep = [keep  keep_partition];
end

numkeep = length(keep); %to adjust for rounding off errors not balancing out
keep = sort(keep);

thetaskeep = thetas(keep);
Pgiven = P(:, keep);

distmatrix = zeros(numkeep);
for i = 1:numkeep
    for j = 1:numkeep
        distmatrix(i,j) = mydist(Pgiven(:,i), Pgiven(:,j));
    end
end
dlmwrite('distances.csv', [thetaskeep; distmatrix]);

expectedanswer = (thetaskeep - thetaskeep(1))';
expectedanswer  = expectedanswer * mode(sign(expectedanswer));

%%%% Add noise

ref = std(Pgiven(:));
sigmaNoise = sigmaNoiseFraction * ref;
noise = normrnd(0, sigmaNoise, size(Pgiven));
Pgiven = Pgiven + noise;
Pgiven = max(0, Pgiven);


noisyoriginal_distmatrix = zeros(numkeep);
for i = 1:numkeep
    for j = 1:numkeep
        noisyoriginal_distmatrix(i,j) = mydist(Pgiven(:,i), Pgiven(:,j));
    end
end


%%%% Calculate Projection moments from projections %%%%

kmax = numkeep - 1; %Number of projection moments. This equality is because this kmax is the number needed for calculating image moments
%"The image moments of the kth order are determined by 
% a set of projection moments of the kth order from
% k+1 given views."

%TODO: Only first order moments are needed
M=ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
for k = 0:kmax
    for i = 1:size(Pgiven,2)
        M(i, k+1) = calculateProjectionMoment(Pgiven(:,i), svector, k);
    end
end

%The resulting M is the limited set of Projection Moments
% M(thetanumber, momentnumber)
      M(:,1); % 0th projection moment
      M(:,2); %

%%%% Calculate Real image moments
IMbasetruth = [ imageMomentFromImage(I, 1, 0, size(I)/2, smax);
                imageMomentFromImage(I, 0, 1, size(I)/2, smax); 
              ];
      
      
%% %%%%%%%%%%% B. (POTENTIALLY) DENOISING TO BE DONE HERE %%%%%%%%

      
%% %%%%%%%%%%%%%%% C. RECOVER ANGLES %%%%%%%%%%%%%%%%%
M_original = M;
expectedanswer_original = expectedanswer;
    %These will both be reordered

     
candidatebesterror = zeros(numkeep - 1, 1);
IMvectors = zeros(numkeep, 2,3); %index1 = candidate_number
                                %index2 = (1,0) | (0,1)
                                %index3 = basetruth | assumed | reconstructed
                                
thetaspredicted = zeros(numkeep, numkeep); %index1: theta2 candidate. Index2: theta_i
thetaerrorvector = zeros(numkeep, 1);


for candidate_index = 2:numkeep %candidate for theta2
    disp(strcat('Considering angle', int2str(candidate_index), ' for theta2'));
    
    expectedanswer = expectedanswer_original;
    expectedanswer(2) = expectedanswer_original(candidate_index);
    expectedanswer(candidate_index) = expectedanswer_original(2);
    
    M = M_original;
    M(2, :) = M_original(candidate_index, :);
    M(candidate_index, :) = M_original(2, :);
    
    %%% Full procedure for every value for theta2
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
    
    
    IMvectors(candidate_index, :, 1) = IMbasetruth;
    IMvectors(candidate_index, :, 2) = bestIMassumed;
    IMvectors(candidate_index, :, 3) = bestIMreconstructed;
    
            %index1 = candidate_number
            %index2 = (1,0) | (0,1)
            %index3 = basetruth | assumed | reconstructed
    thetaspredicted(candidate_index, :) = bestanswer;
    temp = thetaspredicted(candidate_index, 2);                                  %reversing the swap
    thetaspredicted(candidate_index, 2) = thetaspredicted(candidate_index, candidate_index);      %reversing the swap
    thetaspredicted(candidate_index, candidate_index) = temp;                    %reversing the swap
    
    thetaerrorvector(candidate_index) = min( sum(abs(bestanswer - expectedanswer)), sum(abs(bestanswer + expectedanswer))) ;

end

% disp(strcat('Expected Answer | BestAnswer')); 
% [expectedanswer, bestanswer] %#ok<NOPTS>
%disp(strcat('Total Absolute Error: ',num2str(sum(abs(bestanswer - expectedanswer)))));

% figure();
% plot(rangethetas, totalerrorvector);

% low  = max(1, bestanswer(2) - 30);
% high = min(bestanswer(2) + 30, 179);
% figure();
% plot(rangethetas(low:high), totalerrorvector(low:high));
%plot(rangethetas(1:150), totalerrorvector(1:150));

%% %%%%%%%%%%%%%%% D. IDENTIFYING BAD THETAS %%%%%%%%%%%%%%%%%

%%Direct median
% bestthetaspredicted_median = (median(thetaspredicted))';
% bestthetaspredicted_median = bestthetaspredicted_median * mode(sign(bestthetaspredicted_median));
% %[expectedanswer_original, bestthetaspredicted_median, abs(bestthetaspredicted_median - expectedanswer_original) <= 5 ]
% 
% [bestthetaspredictedsorted_median, order_median] = sort(bestthetaspredicted_median);
%  
% expectedanswer_original_ordered = expectedanswer_original(order_median);
% mediancorrect = abs(bestthetaspredictedsorted_median - expectedanswer_original_ordered) <= 5;
% 
% format shortG
% disp('Actual | median | mediancorrect?');
% [expectedanswer_original_ordered, bestthetaspredictedsorted_median, mediancorrect]
% sum(1 - (abs(bestthetaspredictedsorted_median - expectedanswer_original_ordered) <=5))


%% convert to 0 to 180, then take median

thetaspredicted = mod(thetaspredicted, 180);

bestthetaspredicted_median = (median(thetaspredicted))';
bestthetaspredicted_median = bestthetaspredicted_median * mode(sign(bestthetaspredicted_median));
%[expectedanswer_original, bestthetaspredicted_median, abs(bestthetaspredicted_median - expectedanswer_original) <= 5 ]

[bestthetaspredicted_median_sorted, order_median] = sort(bestthetaspredicted_median);
 
expectedanswer_original_reordered = expectedanswer_original(order_median);
mediancorrect = abs(bestthetaspredicted_median_sorted - expectedanswer_original_reordered) <= 5;

format shortG
disp('Actual | median | mediancorrect?');
[expectedanswer_original_reordered, bestthetaspredicted_median_sorted, mediancorrect]
sum(1 - (abs(bestthetaspredicted_median_sorted - expectedanswer_original_reordered) <=5))



%% Decide good and bad thetas as per inter-angle distances
%Order as per tsp
[tsporder1, tsporder2] = orderbytsp(Pgiven);

%Order as per deduced angles
deducedorder = order_median;

if sum(abs(deducedorder - tsporder1)) < sum(abs(deducedorder - tsporder2))
    tsporder = tsporder1;
else
    tsporder = tsporder2;
end

[deducedorder, tsporder, abs(deducedorder - tsporder)]

ordercorrectness = ones(numkeep, 1);
for i = 1:numkeep
    if abs(i - find(deducedorder(i) == tsporder)) > 2
        ordercorrectness(i) = 0;
    end
end

ordercorrectness = ones(numkeep, 1);
defecit = 0;
for i = 1:numkeep
    diff = i - find(deducedorder(i) == tsporder);
    if diff == defecit
        continue;
    end
    if abs(diff - defecit) == 1
        defecit = diff;
        continue
    end
    if abs(diff - defecit) > 1
        ordercorrectness(i) = 0;
        if diff > 0
            % actual angle is further down
            defecit = defecit + 1; %adding 1 to a negtive number
        else
            % actual angle should have been up in the list
           defecit = defecit + 0; % do nothing. making it explicit.
        end
        
    end
end


[deducedorder, tsporder, abs(tsporder- deducedorder), ordercorrectness]

[bestthetaspredicted_median_sorted, order_median] = sort(bestthetaspredicted_median);
expectedanswer_original_reordered = expectedanswer_original(order_median);
mediancorrect = abs(bestthetaspredicted_median_sorted - expectedanswer_original_reordered) <= 3;
% disp('Predicted | Actual | Do we think it is correct? | Is it actually correct?');
% [bestthetaspredicted_median_sorted, expectedanswer_original_reordered, ordercorrectness, mediancorrect]

disp('Predicted | Actual | Correct?');
[bestthetaspredicted_median_sorted, expectedanswer_original_reordered, mediancorrect]


%% %%%%%%%%%%%%%%%%% IMAGE RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%

% Using ORIGINAL angles
% % I_reconstructed = iradon(Pgiven, thetaskeep);
% % figure();
% % imshow(I_reconstructed,'InitialMagnification','fit');


% Using ALL RECONSTRUCTED angles (including incorrect ones)
% % Pgiven_reordered = Pgiven(order_median)


% a(logical([1,0,1,1,0,0,0,1,0,0])) = []
