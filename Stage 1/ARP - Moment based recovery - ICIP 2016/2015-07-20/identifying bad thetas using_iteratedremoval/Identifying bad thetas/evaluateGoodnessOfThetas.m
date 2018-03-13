function [IMreconstructed, error] = evaluateGoodnessOfThetas(estimatedThetas, IMassumed, allPM)
    % Given a set of thetas, estimates how good the fit is by comparing
    % reconstructed IM with original IM. 
    % IM is generated using a least squares approach
    %absolute error is used as it seems to be performing much better
    %(robust to noise) than squared error
    
    A = [cosd(estimatedThetas), sind(estimatedThetas)];
    IMreconstructed = A \ allPM;
    % A * IM = PM
    
    error = sum(abs(IMassumed - IMreconstructed));
    %error = sum((IMassumed - IMreconstructed) ) .^ 2);
    
end
    
