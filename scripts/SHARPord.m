function [refinedProjections, thetasestimated, shiftsestimated] = ...
    SHARPord(Pgiven, svector, sigmaNoise, shift_amplitude)
    numkeep = size(Pgiven, 2);
    kmax = numkeep - 1;
    numstarts = 15;
    Ord = 8;

    Pgiven = denoise(Pgiven, sigmaNoise, 50, 700);
    Pgiven = max(0, Pgiven);
    refinedProjections = Pgiven;

    numiter = 30;
    Ehlccvalues = zeros(numiter * (numkeep - 1), 1);   
    deltas = zeros(numiter, 1);

    thetasestimated_bystart = zeros(numkeep, numstarts);
    shiftsestimated_bystart = zeros(numkeep, numstarts);
    Ehlccvalues_bystart = zeros(numstarts, 1);

    for start = 1:numstarts
        % Initialize shifts randomly.
        shiftestimated = randi([-shift_amplitude shift_amplitude], numkeep, 1);
        shiftedPgiven = zeros(size(Pgiven));
        for i=1:size(Pgiven, 2)
            shiftedPgiven(:, i) = circshift(Pgiven(:,i), shiftestimated(i));
        end
        PMord = ...
            assemblePMord(shiftedPgiven, kmax, svector, Ord, numkeep);
        
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
                finalbesttheta = thetasestimated(i);
                bestshift = shiftestimated(i);
                bestshiftE = Ehlcc;
                bestE = Ehlcc;
                
                for s = -shift_amplitude:shift_amplitude
                    shift_iter = shiftestimated;
                    shift_iter(i) = s;
                    shiftedPgiven = zeros(size(Pgiven));
                    for k=1:size(Pgiven, 2)
                        shiftedPgiven(:, k) = circshift(Pgiven(:,k), shift_iter(k));
                    end
                    PMord = ...
                        assemblePMord(shiftedPgiven, kmax, svector, Ord, numkeep);
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
                    
                    if bestE < bestshiftE
                        bestshift = s;
                        bestshiftE = bestE;
                        finalbesttheta = besttheta;
                    end
                end

                thetasestimated(i) = finalbesttheta;
                shiftestimated(i) = bestshift;
                
                shiftedPgiven = zeros(size(Pgiven));
                for k=1:size(Pgiven, 2)
                    shiftedPgiven(:, k) = circshift(Pgiven(:,k), shiftestimated(k));
                end
                PMord = ...
                    assemblePMord(shiftedPgiven, kmax, svector, Ord, numkeep);

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
        shiftsestimated_bystart(:, start) = shiftestimated;
    end % multistart loop

    [~, optstart] = min(Ehlccvalues_bystart);
    thetasestimated = thetasestimated_bystart(:, optstart);
    shiftsestimated = shiftsestimated_bystart(:, optstart);
end
