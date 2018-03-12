function thetasestimated = ...
    moment_angle_estimation(projections, thetasestimated)
    I_current = iradon(projections, thetasestimated);
    I_current = I_current(2:end-1, 2:end-1); % crop

    E_current = (norm(projections - radon(I_current, thetasestimated),2)) ^ 2;
    E_previter = E_current;

    numiter = 30;
    E_values = zeros(numiter, 1);

    ctr = 0;
    for iteration = 1:numiter
        ctr = ctr + 1;
        for i = 1:size(thetasestimated, 2)
            besttheta = thetasestimated(i);
            bestE = E_current;

            for t = -0:0.1:180
                thetas_iter = thetasestimated;
                thetas_iter(i) = t;

                E_t = (norm(projections - radon(I_current, thetas_iter),2)) ^ 2;

                if E_t < bestE
                    besttheta = t;
                    bestE = E_t;
                end
            end
            thetasestimated(i) = besttheta;
            I_current = iradon(projections, thetasestimated);
            I_current = I_current(2:end-1, 2:end-1); % crop
            E_current = (norm(projections - radon(I_current, thetasestimated),2)) ^ 2;
        end
        E_iter = E_current;
        E_values(ctr) = E_iter;
        delta = E_previter - E_iter;
        if delta < 0.0001
            break
        end
        E_previter = E_iter;
    end
end
