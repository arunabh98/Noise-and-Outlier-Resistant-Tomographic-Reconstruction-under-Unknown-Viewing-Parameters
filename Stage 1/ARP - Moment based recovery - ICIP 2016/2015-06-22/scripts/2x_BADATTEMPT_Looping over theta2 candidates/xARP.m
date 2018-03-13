%%
% Reconstruction from projection data from unknown angles
% In this script we'll implement angle recovery. From there, reconstruction
% is quite directly possible

clear; clc;
%rng(50);

%% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
thetas = (0:179); % Projections to sample from
numkeep = 10; % Number of angles to recompute from
sigmaNoiseFraction = 0.0005;


%% %%%%%%%%%%%%% A. GENERATE INPUT DATA %%%%%%%%%%%%%%%%%

%%%% Define image I %%%%
I = phantom(200);
%I = rgb2gray(imread('../images/640px-napoleon.jpg'));
 
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
keep = randperm(size(P,2), numkeep);
%keep = randi([1 size(P,2)], numkeep, 1);
%keep = sort(keep);


thetaskeep = thetas(keep);
Pgiven = P(:, keep);
expectedanswer = (thetaskeep - thetaskeep(1))';

%%%% Add noise

ref = std(Pgiven(:));
sigmaNoise = sigmaNoiseFraction * ref;
Pgiven = Pgiven + normrnd(0, sigmaNoise, size(Pgiven));
Pgiven = max(0, Pgiven);


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
          
%% %%%%%%%%%%% (POTENTIALLY) DENOISING TO BE DONE HERE %%%%%%%%

      
%% %%%%%%%%%%%%%%% B. RECOVER ANGLES %%%%%%%%%%%%%%%%%
Pgiven_original = Pgiven;
expectedanswer_original = expectedanswer;
    %These will both be reordered
    
candidatebesterror = zeros(numkeep - 1, 1);

IMvectors = zeros(numkeep, 2,3); %index1 = candidate_number
                                %index2 = (1,0) | (0,1)
                                %index3 = basetruth | assumed | reconstructed

for candidate_index = 2:5%numkeep %candidate for theta2
    disp(strcat('Considering angle', int2str(candidate_index), ' for theta2'));
    expectedanswer = expectedanswer_original;
    expectedanswer(2) = expectedanswer_original(candidate_index);
    expectedanswer(candidate_index) = expectedanswer_original(2);
    
    Pgiven = Pgiven_original;
    Pgiven(2) = Pgiven_original(candidate_index);
    Pgiven(candidate_index) = Pgiven_original(2);
    
    IMvectors(candidate_index, :, 1) = IMbasetruth;
    
    %%% Full procedure for every value for theta2
    theta1 = 0;
    bestanswer = zeros(numkeep);
    besterror = Inf;

    rangethetas = 1:179;
    ctr = 1;
    totalerrorvector = zeros(size(rangethetas));

    %% Rearrange thetas to put appropriate angle in position (2)

    for theta2 = rangethetas    %TODO: Increase resolution
                                %     (Possibly in multiple iterations
                                %      in negihbourhood of best estimate)
        
        %%  Calculate (first order) image moments from 2 projection moments and angles

        PM1 = M(:, 2); % 1 indicates order 

        A = [cosd(theta1), sind(theta1); cosd(theta2), sind(theta2)];
        IMassumed = A \ PM1(1:2);   % = A^-1 * PM, using first two values
        IMvectors(candidate_index,:,2) = IMassumed;
                                   %2 => assumed

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
            
            IMvectors(candidate_index,:,3) = IMreconstructed;
            IMassumed
        end

    end
    
    
    
end

IMvectors

    disp(strcat('Expected Answer | BestAnswer')); 
    [expectedanswer, bestanswer]
    TAE  = min( sum(abs(bestanswer - expectedanswer)), sum(abs(bestanswer + expectedanswer))) ;
    disp(strcat('Total Absolute Error = ',num2str(TAE)));

    figure();
    plot(rangethetas, totalerrorvector);

    low  = max(1, bestanswer(2) - 30);
    high = min(bestanswer(2) + 30, 179);
    figure();
    plot(rangethetas(low:high), totalerrorvector(low:high));
    plot(rangethetas(1:100), totalerrorvector(1:100));
    

figure();
theta2candidates = thetaskeep(2:numkeep);
%M1,0
plot(theta2candidates, IMvectors(2:numkeep,1,1), 'r*-', ...
     theta2candidates, IMvectors(2:numkeep,1,2), 'go-', ...
     theta2candidates, IMvectors(2:numkeep,1,3), 'b+-');

%M0,1
plot(theta2candidates, IMvectors(2:numkeep,2,1), 'r*', ...
     theta2candidates, IMvectors(2:numkeep,2,2), 'go', ...
     theta2candidates, IMvectors(2:numkeep,2,3), 'b+');

 
 
%% %%%%%% COMMENTS %%%%%%%%
%1. Looping over 360 degrees becomes necessary because theta2 - theta1
%     could be as large as 180 deg
%2. Guessing and fixing one theta allows first order moments. Two loops
%     will allow 2nd order moments. We can trade-off between processing
%     time and accuracy. How is overall process accuracy affected by this?
% 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%  disp('Image moments directly from image')
%  imageMomentFromImage(I, 1, 0, [100;100], smax)
%  imageMomentFromImage(I, 0, 1, [100;100], smax)








