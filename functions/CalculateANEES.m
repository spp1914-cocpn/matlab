function anees = CalculateANEES(trueStates, estimates, estimateCovs)

    s1 = size(trueStates); % first dim: state, second dim: numTimeSteps, --third dim: numRuns
    
    if isequal(s1, size(estimates))
        numTimeSteps = s1(2);
        %if numel(s1) == 3
        %    numRuns = s1(3);
        %else
        %    numRuns = 1;
        %end
        anees = zeros(numTimeSteps, 1);
        errors = estimates - trueStates;
        for k=1:numTimeSteps
            nees = 0;
            %for r=1:numRuns
                error = errors(:, k);
                nees = nees + (error' / estimateCovs(:, :, k)) * error;
            %end
            anees(k) = nees;
        end
        anees = anees / (s1(1));
    else
        anees = [];
    end
end

