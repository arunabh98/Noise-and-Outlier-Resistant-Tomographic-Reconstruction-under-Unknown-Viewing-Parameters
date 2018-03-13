



distm = normrnd(0,10,10);
for i = 1:10
    for j = 1:10
        distm(i,j) = abs(distm(i,j));
        distm(j,i) = distm(i,j);
    end
    distm(i,i) = 0
end

xy = abs(normrnd(0,10,10,5));

%userConfig = struct('popSize',200,'numIter',1e4); 



