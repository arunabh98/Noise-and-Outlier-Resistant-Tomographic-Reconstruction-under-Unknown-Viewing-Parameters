function [order1, order2, dist] = orderbytsp(Pmatrix)
    
    userConfig = struct('xy', Pmatrix', 'showProg',false,'showResult',false,'showWaitbar',true);
    resultStruct = tsp_ga(userConfig);
    
    %Reorder to start at <1>
    order1 = resultStruct.optRoute;

    start = find(order1 == 1);
    if start ~= 1
        order1 = [order1(start:length(order1)), order1(1:(start-1))];
    end
    order2 = [1, flip(order1(2:length(order1)))];
    %isequal(order1, 1:numkeep) | isequal(order2, 1:numkeep)
    order1 = order1';
    order2 = order2';
    dist = resultStruct.minDist;
    
end


