clear; clc;
%%%% Define image I %%%%
%rng(22);

% Ibeam
I = zeros(100,100);
x1 = 20; y1 = 10;
x2 = 80; y2 = 30;
x3 = 40; y3 = 30;
x4 = 60; y4 = 70;
x5 = 10; y5 = 70;
x6 = 90; y6 = 90;

I(y1:y2, x1:x2) = 1;
I(y3:y4, x3:x4) = 1;
I(y5:y6, x5:x6) = 1;

%Display image
imshow(I,'InitialMagnification','fit');

%%%% Calculate projections %%%%
%%%Specific theta % For testing only
    theta = 0;
    Ptheta = radon(I, theta);
    figure();
    plot(Ptheta);
    xlabel('s');
    ylabel('P_0(s)');
    title('Projections at \theta = 0^o');

    theta = 90;
    Ptheta = radon(I, theta);
    figure();
    plot(Ptheta);
    xlabel('s');
    ylabel('P_{90}(s)');
    title('Projections at \theta = 90^o');

%All thetas:

thetas = (0:179);
[P, svector] = radon(I, thetas); % P(s, theta)
figure();
imshow(P / max(max(P)), [], 'Xdata', [0 179] , 'Ydata', svector, 'InitialMagnification','fit');
iptsetpref('ImshowAxesVisible','on')
xlabel('Angle \theta')
ylabel('s')
title('Projections P_\theta(s)')

% Normalize s to a unit circle
smax = max(abs(svector));
svector = svector / smax;
%

gapsvector = [170 85 34 17 10 5 2];

IMerrors = zeros(length(gapsvector), 1);
PMerrors = zeros(length(gapsvector), 1);

%ctr = 1
%for gap = gapsvector
gap = 10
    
%%%% Create input data %%%%
numkeep = 12;
keep = randperm(size(P,2), numkeep);
%keep = randi([1 size(P,2)], numkeep, 1);
keep = sort(keep);
Pgiven = P(:, keep);
thetasgiven = thetas(keep);

kmax = numkeep - 1; %Number of projection moments. This equality is because this kmax is the number needed for calculating image moments
%"The image moments of the kth order are determined by 
% a set of projection moments of the kth order from
% k+1 given views."

M=ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
for k = 0:kmax
    for i = 1:size(Pgiven,2)
        M(i, k+1) = calculateProjectionMoment(Pgiven(:,i), svector, k);
    end
end

%size(M)
%The resulting M is the limited set of Projection Moments
% M(thetanumber, momentnumber)
%       M(:,1);
%      M(:,2)
%      M(:,10)
%      M(:,14)
% 
%      M


%%%% Calculate Image moments %%%%

IM = zeros(kmax+1, kmax+1);
%This stores all image moments, column wise. However, since there are
%only k+1 image moments of order k, the rest of the rows are zero
%IM(p-number, momentnumber)

for k = 0:kmax
    %PMk: Projection moments of order k
    %IMk: Image moments of order k
    %Note:  PMk is a vector because there are multiple views involved
    %       IMk is a vector because IMk has k+1 values, one for each l

    PMk = M(1:(k+1),k+1);
    A = assembleA(thetasgiven, k); %  Assembling the A matrix [ IM(k) = A^-1 PM(k) ]
    inv(A);
    
    IMk = A \ PMk; %This is basically IMk = inv(A) * PMk; %If incidentally inv(A) is not defined, add noise
    IM(1:(k+1),k+1) = IMk;
end

%%%% Verification of IMs
IMactual = zeros(k+1,k+1);
center = floor((size(I)+1)/2);

for k = 0:kmax
    for p = k:-1:0
        IMDpq = imageMomentFromImage(I, p, k-p, center, smax);
        IMactual(k-p+1,k+1)= IMDpq;
    end
end
dlmwrite('output/ImageMomentsConstructed.csv', IM);
dlmwrite('output/ImageMomentsActual.csv', IMactual);

IMerror = abs(IM - IMactual) %Cannot take relative errors because IMactual has values very very close to 0 (sometimes 0)
               
       
%%%%%%%%%%%

%%%% Produce Projection moment of arbitrary (new) thetas  %%%%
% We can verify that correct PMs are being computed at this point. 
%(We have the original Ps stored! So we can get the PMs)

newthetas = 0:gap:170;  %randomly, say   

PMreconstructed = zeros(length(newthetas), kmax+1);
for k = 0:kmax
    for t = 1:length(newthetas)
        mockthetas = zeros(k+1);
        mockthetas(1) = newthetas(t);
        A = assembleA(mockthetas, k);
        A = A(1,:);
        %Now, PM = A * IM
        PMnew = A * IM(1:k+1, k+1);
        PMreconstructed(t, k+1) = PMnew;
    end
end

PMactual = zeros(length(newthetas), kmax+1);
for k = 0:kmax
    for t = 1:length(newthetas)
        PMktheta = calculateProjectionMoment(P(:,newthetas(t)+1), svector, k);
        PMactual(t, k+1) = PMktheta;
    end
end
PMerror = abs(PMreconstructed - PMactual)
%PMerrors(ctr) = median(median(PMerror));
%ctr = ctr + 1
% end
% plot(170 ./ gapsvector + 1, PMerrors)
% xlabel('#Reconstructed views');
% ylabel('Median error');
% title('Median Reconstruction error ariation with number of reconstructed views');



dlmwrite('output/ProjectionMomentsConstructed.csv', PMreconstructed);
dlmwrite('output/ProjectionMomentsActual.csv', PMactual);

densityvector = svector; %-1:0.005:1;


P_reconstructed = zeros(length(densityvector), length(newthetas));


for t = 1:length(newthetas)
    Ptheta = zeros(length(svector), 1);
    theta = newthetas(t);
    for i = 1:length(densityvector)
        s = densityvector(i);
        Pstheta = reconstructPfromM(s, kmax, PMreconstructed(t, :));
        Ptheta(i) = Pstheta;
    end
    P_reconstructed(:,t) = Ptheta;
    %figure();
    %plot(Ptheta)
end

P_reconstructed2 = P_reconstructed / max(max(P_reconstructed));
Ireconstructed = iradon(P_reconstructed, newthetas);
Ireconstructed = iradon(P_reconstructed, newthetas);
Ireconstructed = iradon(Pgiven, thetasgiven);
figure();
imshow(Ireconstructed,'InitialMagnification','fit');










    
    