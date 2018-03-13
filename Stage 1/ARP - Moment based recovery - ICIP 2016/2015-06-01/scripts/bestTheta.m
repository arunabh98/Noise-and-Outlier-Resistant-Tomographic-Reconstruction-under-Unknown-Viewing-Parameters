function [bestt, besterror] = bestTheta(IM, PM, theta1)
    if ~exist('theta1', 'var')
        theta1 = 0;
    end
        
    bestt = 0;
    besterror = Inf;
    for t = [-179:-1, 1:179]
        A = [cosd(theta1), sind(theta1); cosd(t), sind(t)];
        IMcalculated = A \ PM;
        error = sum(sum(abs(IMcalculated - IM)));
        
        if error < besterror
            bestt = t;
            besterror = error;
        end
    end
        
end
