%%
% Reconstruction from projection data from unknown angles
% In this script we'll implement angle recovery. From there, reconstruction
% is quite directly possible

clear; clc;
rng(100);

%% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
thetas = (0:179); % Projections to sample from
numkeep = 30; % Number of angles to recompute from
sigmaNoiseFraction = 0.001;


%% %%%%%%%%%%%%% A. GENERATE INPUT DATA %%%%%%%%%%%%%%%%%

%%%% Define image I %%%%
I = phantom(200);
%I = rgb2gray(imread('../images/640px-napoleon.jpg'));
%I = rgb2gray(imread('../images/256px-Earthrise.jpg'));
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
keep = sort(keep);


thetaskeep = thetas(keep);
Pgiven = P(:, keep);
expectedanswer = (thetaskeep - thetaskeep(1))';

%%%% Add noise

ref = std(Pgiven(:));
sigmaNoise = sigmaNoiseFraction * ref;
noise = normrnd(0, sigmaNoise, size(Pgiven));
Pgiven = Pgiven + noise;
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

bestthetaspredicted = (median(thetaspredicted))';

% Take one theta at a time.  Construct Im from just that one, vs from all
% the others. If this difference is large, this is a bad theta

%The idea is that if most of the thetas are correct, the least squares
%estimate will be correct


fusionerrorvector = zeros(numkeep, 2);
imlsvector = zeros(numkeep, 2);
imivector  = zeros(numkeep, 2);

for candidate_index = 2:numkeep
    theta_i = thetaskeep(candidate_index);
    otherthetas = thetaskeep';
    otherthetas([1,candidate_index]) = [];
    PMtheta_i = M_original(candidate_index, 2);
    PM0 = M_original(1, 2);
    PMotherthetas = M_original(:, 2);
    PMotherthetas([1,candidate_index]) = [];
    
    
    %[imls, imi] = 
    error = errorInMomentsByTheta(theta_i, otherthetas, PMtheta_i, PM0, PMotherthetas);
    fusionerrorvector(candidate_index, :) = error;
    %imlsvector(candidate_index, :) = imls;
    %imivector(candidate_index, :)  = imi;
end

[expectedanswer_original, bestthetaspredicted, abs(bestthetaspredicted - expectedanswer_original) <= 5, fusionerrorvector , (sum(abs(noise)))']
corr((sum(abs(noise)))', abs(fusionerrorvector(:,1)))
corr((sum(abs(noise)))', abs(fusionerrorvector(:,2)))



%% %%%%%%%%%%%%%%% E. RECONSTRUCTION ANALYSIS %%%%%%%%%%%%%%%%%

IMbasetruth
bestIMassumed
bestIMreconstructed

TAEangles   = min( sum(abs(bestanswer - expectedanswer)), sum(abs(bestanswer + expectedanswer))) ;
disp(strcat(['Total Absolute Error in angles = ', num2str(TAEangles)]));

% apparentError = sum(abs(bestIMreconstructed - bestIMassumed)); %TODO: This is using last iteration
% actualError   = sum(abs(bestIMreconstructed - IMbasetruth));   %TODO: This is using last iteration
% disp(strcat(['Total Absolute Error in Image moments (from assumed)    = ', num2str(apparentError)]));
% disp(strcat(['Total Absolute Error in Image moments (from base truth) = ', num2str(actualError)]));


%% Base truth vs Calculate vs Reconstructed
theta2candidates = expectedanswer_original(2:numkeep);

%figure();
%M1,0
plot(theta2candidates, IMvectors(2:numkeep,1,1), 'r-', ...
     theta2candidates, IMvectors(2:numkeep,1,2), 'go-', ...
     theta2candidates, IMvectors(2:numkeep,1,3), 'b+-');
xlabel('(Real) Value of theta2 candidate')
title('M_{1,0}');
legend('Base truth', 'Calculated', 'Reconstructed', 'Location', 'best');

%figure();
%M0,1
plot(theta2candidates, IMvectors(2:numkeep,2,1), 'r-', ...
     theta2candidates, IMvectors(2:numkeep,2,2), 'go-', ...
     theta2candidates, IMvectors(2:numkeep,2,3), 'b+-');
xlabel('(Real) Value of theta2 candidate')
title('M_{0,1}');
legend('Base truth', 'Calculated', 'Reconstructed', 'Location', 'best');

%% %%% Apparent errors AND actual error

%M1,0
M10apparenterror = abs(IMvectors(2:numkeep,1,3) - IMvectors(2:numkeep,1,2));
M10actualerror   = abs(IMvectors(2:numkeep,1,3) - IMvectors(2:numkeep,1,1));
correlation10 = corrcoef(M10apparenterror, M10actualerror);

%figure();
plot(theta2candidates, M10apparenterror, 'go-', ...
     theta2candidates, M10actualerror,   'ro-');
disp('M_1,0');

correlation10 = num2str(correlation10(1,2));
disp(strcat('Correlation in apparent and actual error: ', correlation10));
title(strcat(['M_{1,0}', '  Correlation: ', correlation10]));
legend('Apparent Error', 'Actual Error', 'Location', 'best');



%M0,1
M01apparenterror = abs(IMvectors(2:numkeep,2,3) - IMvectors(2:numkeep,2,2));
M01actualerror   = abs(IMvectors(2:numkeep,2,3) - IMvectors(2:numkeep,2,1));
correlation01 = corrcoef(M01apparenterror, M01actualerror);
%figure();
plot(theta2candidates, M01apparenterror, 'go-', ...
     theta2candidates, M01actualerror,   'ro-');
disp('M_0,1');

correlation01 = num2str(correlation01(1,2));
disp(strcat('Correlation in apparent and actual error: ', num2str(correlation01)));

legend('Apparent Error', 'Actual Error', 'Location', 'best');
title(strcat(['M_{0,1}', '  Correlation: ', correlation01]))




%% Apparent error AGAINST actual error
%M1,0
[M10apparenterrorsorted, order] = sort(M10apparenterror);
M10actualerrorsorted = M10actualerror(order);
 
%figure();
plot(M10apparenterrorsorted, M10actualerrorsorted, 'ro-');
title(strcat(['M_{1,0} - Apparent vs actual error |', '  Correlation: ', correlation10]));
xlabel('Apparent Error');
ylabel('Actual Error');


%M0,1
[M01apparenterrorsorted, order] = sort(M01apparenterror);
M01actualerrorsorted = M01actualerror(order);
 
%figure();
plot(M01apparenterrorsorted, M01actualerrorsorted, 'ro-');
title(strcat(['M_{1,0} - Apparent vs actual error |', '  Correlation: ', correlation01]));
xlabel('Apparent Error');
ylabel('Actual Error');

%% %%%%%%%%%%%%%%% F. BEST RECONSTRUCTION %%%%%%%%%%%%%%%%%

%% Using M1,0
[besttheta2candidate_error, besttheta2candidate_index]  = min(M10apparenterror); %#ok<ASGLU>
bestthata2candidate1 = expectedanswer_original(besttheta2candidate_index)
% Using M0,1
[besttheta2candidate_error, besttheta2candidate_index]  = min(M01apparenterror);
bestthata2candidate2 = expectedanswer_original(besttheta2candidate_index)


%% %%%%%%%%%%%%%%% G. WRITING OUT %%%%%%%%%%%%%%%%%
% format shortG
% dlmwrite('data.csv', [thetaskeep', thetaerrorvector, (sum(noise))', IMvectors(:,:,1), IMvectors(:,:,2), IMvectors(:,:,3)])
 dlmwrite('thetasdata.csv', [thetaskeep', thetaspredicted']);
% format

%% Reevaluate for selected candidate theta2
% Theoretically, this sould already have been stored, but this is simpler
% for now. TODO: Implement the storing procedure


%% %%%%%% COMMENTS %%%%%%%%
%1. Looping over 360 degrees becomes necessary because theta2 - theta1
%     could be as large as 180 deg
%2. Guessing and fixing one theta allows first order moments. Two loops
%     will allow 2nd order moments. We can trade-off between processing
%     time and accuracy. How is overall process accuracy affected by this?
% 




