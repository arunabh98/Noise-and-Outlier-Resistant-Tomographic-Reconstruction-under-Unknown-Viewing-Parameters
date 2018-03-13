%%
% Reconstruction from projection data from unknown angles
% In this script we'll implement angle recovery. From there, reconstruction
% is quite directly possible

close all; clear; clc;
rng(6);

tic
%% %%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%
thetas = (0:179); % Projections to sample from
numkeep = 101; % Number of angles to recompute from
sigmaNoiseFraction = 0.1;


%% %%%%%%%%%%%%% A. GENERATE INPUT DATA %%%%%%%%%%%%%%%%%

%%%% Define image I %%%%
%I = phantom(200);
%I = rgb2gray(imread('../images/256px-Earthrise.jpg'));
%I = rgb2gray(imread('../images/200px-beach.jpg'));
I= rgb2gray(imread('../images/200px-woody.jpg'));
%  figure();
  imshow(I,'InitialMagnification','fit');



%%%% Calculate projections %%%%

[P, svector] = radon(I, thetas); % P(s, theta)


    % View Image
    % figure();
    % imshow(P / max(max(P)), [], 'Xdata', [0 179] , 'Ydata', svector, 'InitialMagnification','fit');
    % iptsetpref('ImshowAxesVisible','on')
    % xlabel('Angle \theta')
    % ylabel('s')
    % title('Projections P_\theta(s)')

% Normalize s to a unit circle
smax = max(abs(svector));
svector = svector / smax;

% Sample Uniformly
keep = randperm(size(P,2), numkeep);
%keep = sort(keep);

% OR
% Sample from a skewed distribution
% weights = [0.2 0.3 0.12 0.03 0.35]; 
% partitionsize = length(thetas)/length(weights);
% 
% keep = [];
% for i = 1:length(weights)
%     start = partitionsize * (i - 1);
%     numsample = round(numkeep * weights(i));
%     keep_partition = start + randperm(partitionsize, numsample);
%     keep = [keep  keep_partition];
% end
% 
% numkeep = length(keep); %to adjust for rounding off errors not balancing out

thetaskeep = thetas(keep);
Pgiven = P(:, keep);
Pnonoise = Pgiven;

% distmatrix = zeros(numkeep);
% for i = 1:numkeep
%     for j = 1:numkeep
%         distmatrix(i,j) = mydist(Pgiven(:,i), Pgiven(:,j));
%     end
% end
% dlmwrite('distances.csv', [thetaskeep; distmatrix]);

expectedanswer = (thetaskeep - thetaskeep(1))';
expectedanswer  = expectedanswer * mode(sign(expectedanswer));

%%%% Add noise

ref = std(Pgiven(:));
sigmaNoise = sigmaNoiseFraction * ref;
noise = normrnd(0, sigmaNoise, size(Pgiven));
Pgiven = Pgiven + noise;
%Pgiven = max(0, Pgiven); % Not done here to keep the noise model Gaussian.
%Done after denoising


% noisyoriginal_distmatrix = zeros(numkeep);
% for i = 1:numkeep
%     for j = 1:numkeep
%         noisyoriginal_distmatrix(i,j) = mydist(Pgiven(:,i), Pgiven(:,j));
%     end
% end


%%%% Calculate Projection moments from projections %%%%

kmax = numkeep - 1; %Number of projection moments. This equality is because this kmax is the number needed for calculating image moments
%"The image moments of the kth order are determined by 
% a set of projection moments of the kth order from
% k+1 given views."

M = ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions 
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
      
      
%% %%%%%%%%%%% B. DENOISING %%%%%%%%
disp('Before denoising');
Pnoisy = Pgiven;
norm(Pnonoise - Pgiven, 'fro')
psnr_noisy = 10 * log10( max(max(Pnonoise))^2 / (mean(mean(( Pnonoise - Pnoisy ) .^ 2))) )

P_denoised1 = denoise(Pgiven, sigmaNoise, 50, 200);
%P_denoised2 = denoise2(Pgiven, sigmaNoise, 50, 200);
Pgiven = P_denoised1;

disp('After denoising')
norm(Pnonoise - P_denoised1, 'fro')
%norm(Pnonoise - P_denoised2, 'fro')
figure()
plot(Pnonoise(:,1));
hold on;
plot(P_denoised1(:,1), '-r');
plot(Pnoisy(:,1), '-g');
hold off;

psnr_denoised = 10 * log10( max(max(Pnonoise))^2 / (mean(mean(( Pnonoise - P_denoised1 ) .^ 2))) )


%% %%%%%% B2. Final fixes on P, M %%%%%%%%%%%%%%%
Pgiven = max(0, Pgiven);

% Recalculate M
M=ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
for k = 0:kmax
    for i = 1:size(Pgiven,2)
        M(i, k+1) = calculateProjectionMoment(Pgiven(:,i), svector, k);
    end
end

% figure()
% plot([Pnonoise Pgiven]);

%% %%%%%%%%%%%%%%% C. RECOVER ANGLES %%%%%%%%%%%%%%%%%

% Initialize angles randomly
thetasestimated = randi([1 179], numkeep, 1);
%thetasestimated(1) = 0;
thetasfirstestimate = thetasestimated;



PM1 = M(:, 2); % 2 indicates order

numiter = 30;
Ehlccvalues = zeros(numiter * (numkeep - 1), 1);
Ehlcciter = Inf;

deltas = zeros(numiter, 1);

%% Multi-start

numstarts = 10;
thetasestimated_bystart = zeros(numkeep, numstarts);
Ehlccvalues_bystart = zeros(numstarts, 1);

for start = 1:numstarts
    
    % Initialize angles randomly
    thetasestimated = randi([1 179], numkeep, 1);
    
    % Starting estimate
    A = [cosd(thetasestimated), sind(thetasestimated)];
    IMestimated = A \ PM1;

    Ehlccvec = A * IMestimated - PM1;
    Ehlcc = norm(Ehlccvec, 2);
    ctr = 1;

    Ehlccvalues(ctr) = Ehlcc;


    for iteration = 1:numiter
        disp([num2str(start), ', ', num2str(iteration)]);

        for i = 1:numkeep
            ctr = ctr + 1;

            besttheta = thetasestimated(i);
            bestE = Ehlcc;

            for t = [-179:180]
                A(i, :) = [cosd(t), sind(t)];
                IMestimated = A \ PM1;
                E_tvec = A * IMestimated - PM1;
                E_t = norm(E_tvec, 2);

                if E_t < bestE
                    besttheta = t;
                    bestE = E_t;
                end
            end

            thetasestimated(i) = besttheta;


            A = [cosd(thetasestimated), sind(thetasestimated)];
            IMestimated = A \ PM1;

            Ehlccvec = A * IMestimated - PM1;
            Ehlcc = norm(Ehlccvec, 2);
            Ehlccvalues(ctr) = Ehlcc;

        end

        ctr
        Ehlcciter = Ehlccvalues(ctr)
        Ehlccpreviter = Ehlccvalues(ctr - (numkeep)); % -1 iff theta0 is not updated

        delta = Ehlccpreviter - Ehlcciter 
        deltas(iteration) = delta;
        if delta < 0.0001
            break
        end

    end
    Ehlccvalues_bystart(start) = Ehlccvalues(ctr);
    thetasestimated_bystart(:, start) = thetasestimated;
end % multistart loop


figure();
plot(Ehlccvalues_bystart, 'or');

[~, optstart] = min(Ehlccvalues_bystart);
%optstart = 7;
thetasestimated = thetasestimated_bystart(:, optstart);

%thetasestimated = thetasestimated_bystart(:, 9);

%thetasestimated = mod(thetasestimated, 360);
[~, sortingorder] = sort(expectedanswer);
disp('Expected | Converged | Starting | (1) -(2) | (1) + (2)');
%[expectedanswer(sortingorder) thetasestimated(sortingorder)
%(mod(expectedanswer(sortingorder) - thetasestimated(sortingorder), 360)) (mod((expectedanswer(sortingorder) + thetasestimated(sortingorder)), 360))]
[expectedanswer thetasestimated (mod(expectedanswer - thetasestimated, 360)) (mod((expectedanswer + thetasestimated), 360))]
toc

    
    %Checker
    mode((mod(expectedanswer - thetasestimated, 360)))
    mode((mod(expectedanswer + thetasestimated, 360)))
    
    %col 3
    sum(abs((mod(expectedanswer - thetasestimated, 360)) - 280) <= 0)
    %col 4
    sum(abs((mod(expectedanswer + thetasestimated, 360)) - 97) <= 0)


Ehlccvalues_bystart'

A = [cosd(expectedanswer), sind(expectedanswer)];
IMestimated = A \ PM1;
Ehlccvec = A * IMestimated - PM1;
Ehlcc_expectedanswer = norm(Ehlccvec, 2)

figure();
plot(Ehlccvalues(1:ctr), '-b');
hold on;
majors = 1:(numkeep - 1):ctr;
plot(majors,Ehlccvalues(majors) , 'ro')
    %Plot deltas
    %figure();
    %plot(deltas(1:iteration), '-r.');


%%%%%
num = 3;
%A = [cosd(thetasestimated([1:(num-1), (num+1):numkeep])), sind(thetasestimated([1:(num-1), (num+1):numkeep]))];
Awithout = [cosd(thetasestimated), sind(thetasestimated)];
Awithout(num,:) = [];
IM = Awithout \ PM1([1:(num-1), (num+1):numkeep]);   % = A^-1 * PM, using first two values
possibleAngles = -179:180;

errors = zeros(size(possibleAngles));
ctr = 0;
for theta = possibleAngles
    ctr = ctr + 1;
    
    Awith = [cosd(thetasestimated), sind(thetasestimated)];
    Awith(num,:) = [cosd(theta), sind(theta)];

    IMwith = Awith \ PM1;
    
    error = norm(IM - IMwith);
    errors(ctr) = error;
    %disp([num2str(theta), '  ', num2str(error)]);
end

figure();
plot(possibleAngles, errors);
[~, minpoint] = min(errors);
disp([num2str(thetasestimated(num)), '  ', num2str(possibleAngles(minpoint))]);



%%%%%



















%% %%%%%%%%%%%%%%% D. IDENTIFYING BAD THETAS %%%%%%%%%%%%%%%%%

% %%Direct median
% % bestthetaspredicted_median = (median(thetaspredicted))';
% % bestthetaspredicted_median = bestthetaspredicted_median * mode(sign(bestthetaspredicted_median));
% % %[expectedanswer_original, bestthetaspredicted_median, abs(bestthetaspredicted_median - expectedanswer_original) <= 5 ]
% % 
% % [bestthetaspredictedsorted_median, order_median] = sort(bestthetaspredicted_median);
% %  
% % expectedanswer_original_ordered = expectedanswer_original(order_median);
% % mediancorrect = abs(bestthetaspredictedsorted_median - expectedanswer_original_ordered) <= 5;
% % 
% % format shortG
% % disp('Actual | median | mediancorrect?');
% % [expectedanswer_original_ordered, bestthetaspredictedsorted_median, mediancorrect]
% % sum(1 - (abs(bestthetaspredictedsorted_median - expectedanswer_original_ordered) <=5))


%% convert to 0 to 180, then take median
% % 
% % thetaspredicted = mod(thetaspredicted, 180);
% % 
% % bestthetaspredicted_median = (mean(thetaspredicted))';
% % bestthetaspredicted_median = bestthetaspredicted_median * mode(sign(bestthetaspredicted_median));
% % %[expectedanswer_original, bestthetaspredicted_median, abs(bestthetaspredicted_median - expectedanswer_original) <= 5 ]
% % 
% % [bestthetaspredicted_median_sorted, order_median] = sort(bestthetaspredicted_median);
% %  
% % expectedanswer_original_reordered = expectedanswer_original(order_median);
% % mediancorrect = abs(bestthetaspredicted_median_sorted - expectedanswer_original_reordered) <= 5;
% % 
% % format shortG
% % disp('Actual | median | mediancorrect?');
% % [expectedanswer_original_reordered, bestthetaspredicted_median_sorted, mediancorrect]
% % sum(1 - (abs(bestthetaspredicted_median_sorted - expectedanswer_original_reordered) <=5))
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % %% Decide good and bad thetas as per inter-angle distances
% % %Order as per tsp
% % [tsporder1, tsporder2] = orderbytsp(Pgiven);
% % 
% % %Order as per deduced angles
% % deducedorder = order_median;
% % 
% % if sum(abs(deducedorder - tsporder1)) < sum(abs(deducedorder - tsporder2))
% %     tsporder = tsporder1;
% % else
% %     tsporder = tsporder2;
% % end
% % 
% % [deducedorder, tsporder, abs(deducedorder - tsporder)]
% % 
% % ordercorrectness = ones(numkeep, 1);
% % for i = 1:numkeep
% %     if abs(i - find(deducedorder(i) == tsporder)) > 2
% %         ordercorrectness(i) = 0;
% %     end
% % end
% % 
% % ordercorrectness = ones(numkeep, 1);
% % defecit = 0;
% % for i = 1:numkeep
% %     diff = i - find(deducedorder(i) == tsporder);
% %     if diff == defecit
% %         continue;
% %     end
% %     if abs(diff - defecit) == 1
% %         defecit = diff;
% %         continue
% %     end
% %     if abs(diff - defecit) > 1
% %         ordercorrectness(i) = 0;
% %         if diff > 0
% %             % actual angle is further down
% %             defecit = defecit + 1; %adding 1 to a negtive number
% %         else
% %             % actual angle should have been up in the list
% %            defecit = defecit + 0; % do nothing. making it explicit.
% %         end
% %         
% %     end
% % end
% % 
% % 
% % [deducedorder, tsporder, abs(tsporder- deducedorder), ordercorrectness]
% % 
% % [bestthetaspredicted_median_sorted, order_median] = sort(bestthetaspredicted_median);
% % expectedanswer_original_reordered = expectedanswer_original(order_median);
% % mediancorrect = abs(bestthetaspredicted_median_sorted - expectedanswer_original_reordered) <= 3;
% % disp('Predicted | Actual | Do we think it is correct? | Is it actually correct?');
% % [bestthetaspredicted_median_sorted, expectedanswer_original_reordered, ordercorrectness, mediancorrect]
% % %[expectedanswer_original, bestthetaspredicted_median, abs(expectedanswer_original - bestthetaspredicted_median) <= 3 , ordercorrectness]
% % 
% % 

%% %%%%%%%%%%%%%%%%% IMAGE RECONSTRUCTION %%%%%%%%%%%%%%%%%%%%%%%%

%Using ORIGINAL angles
I_reconstructed = iradon(Pgiven, thetaskeep);
figure();
imshow(I_reconstructed/max(max(I_reconstructed )), 'InitialMagnification','fit');


%Using ALL RECONSTRUCTED angles (including incorrect ones)
I_reconstructed = iradon(Pgiven, thetasestimated);
figure();
imshow(I_reconstructed/max(max(I_reconstructed )),'InitialMagnification','fit');

% a(logical([1,0,1,1,0,0,0,1,0,0])) = []
