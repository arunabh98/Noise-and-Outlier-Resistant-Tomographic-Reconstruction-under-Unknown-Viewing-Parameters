clearfunction [thetasestimated] = ARPord_Kmeans(Pgiven, svector)
    numkeep = size(Pgiven, 2);
    kmax = numkeep - 1;
    numstarts = 3;
    Ord = 8;

    M=ones(size(Pgiven,2), kmax+1); % each column represents kth moments for one k, for all directions
    for k = 0:kmax
        for i = 1:size(Pgiven,2)
            M(i, k+1) = calculateProjectionMoment(Pgiven(:,i), svector, k);
        end
    end
    PMord = reshape(M(:,2:(Ord+1)), Ord * numkeep, 1);
    numiter = 30;
    Ehlccvalues = zeros(numiter * (numkeep - 1), 1);   
    deltas = zeros(numiter, 1);

    thetasestimated_bystart = zeros(numkeep, numstarts);
    Ehlccvalues_bystart = zeros(numstarts, 1);

    for start = 1:numstarts
        % Initialize angles randomly
        thetasestimated = randi([1 179], numkeep, 1);

        A = assembleA(thetasestimated, Ord); %<>%
        IMestimated = A \ PMord;

        Ehlccvec = A * IMestimated - PMord; 
        Ehlcc = norm(Ehlccvec, 2);
        ctr = 1;

        Ehlccvalues(ctr) = Ehlcc;

        for iteration = 1:numiter
            disp([num2str(start), ', ', num2str(iteration)]);

            for i = 1:numkeep
                ctr = ctr + 1;

                besttheta = thetasestimated(i);
                bestE = Ehlcc;

                for t = -179:180

                    thetas_iter = thetasestimated;
                    thetas_iter(i) = t;
                    A = assembleA(thetas_iter, Ord);
                    IMestimated = A \ PMord;
                    E_tvec = A * IMestimated - PMord;
                    E_t = norm(E_tvec, 2);

                    if E_t < bestE
                        besttheta = t;
                        bestE = E_t;
                    end
                end

                thetasestimated(i) = besttheta;

                A = assembleA(thetasestimated, Ord);
                IMestimated = A \ PMord;

                Ehlccvec = A * IMestimated - PMord;
                Ehlcc = norm(Ehlccvec, 2);
                Ehlccvalues(ctr) = Ehlcc;

            end
            Ehlcciter = Ehlccvalues(ctr);
            Ehlccpreviter = Ehlccvalues(ctr - (numkeep)); % -1 iff theta0 is not updated

            delta = Ehlccpreviter - Ehlcciter;
            deltas(iteration) = delta;
            if delta < 0.0001
                break
            end

        end
        Ehlccvalues_bystart(start) = Ehlccvalues(ctr);
        thetasestimated_bystart(:, start) = thetasestimated;
    end % multistart loop

    [~, optstart] = min(Ehlccvalues_bystart);
    thetasestimated = thetasestimated_bystart(:, optstart);
end
