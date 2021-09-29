function nci = CalculateNCI(trueStates, estimates, estimateCovs)
    % Computes the noncredibility index (NCI) of an estimator
    % based on the results of a Monte-Carlo simulation.
    %
    % Literature: 
    %   Li, X. Rong, Zhanlue Zhao, and Vesselin P. Jilkov, 
    %   Practical measures and test for credibility of an estimator, 
    %   Proc. Workshop on Estimation, Tracking, and Fusion—A Tribute to Yaakov Bar-Shalom, 2001.
    %
    %   Li, X. Rong, Zhanlue Zhao, and Vesselin P. Jilkov, 
    %   Estimator’s credibility and its measures,
    %   Proc. IFAC 15th World Congress. 2002.
    
    s1 = size(trueStates); % first dim: state, second dim: numTimeSteps, third dim: numRuns
    
    if isequal(s1, size(estimates))
        numTimeSteps = s1(2);
        if numel(s1) == 3
            numRuns = s1(3);
        else
            numRuns = 1;
        end
        nci = zeros(numTimeSteps, 1);
        scaleFactor = 10 / numRuns;
        errors = estimates - trueStates;
        for k=1:numTimeSteps
            % first, compute the actual sample MSE matrix per time step
            mseMatrix = zeros(s1(1));
            for r = 1:numRuns
               error = errors(:, k, r);
               mseMatrix = mseMatrix + error * error';
            end
            mseMatrix = mseMatrix / r;
            nci(k) = sum(arrayfun(@(r) log10((errors(:, k, r)' / estimateCovs(:, :, k, r)) *  errors(:, k, r)) ...
                - log10((errors(:, k, r)' / mseMatrix) *  errors(:, k, r)), 1:numRuns)) * scaleFactor;
        end
    else
        nci = [];
    end
end

