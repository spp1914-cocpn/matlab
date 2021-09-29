classdef UnitUpdater < Updater
    methods (Access = public)
        function newTransitionMatrix = updateTransitionProbabilityMatrix(this, ...
                measurement, oldTransitionMatrix)
            newTransitionMatrix = eye( length( oldTransitionMatrix ) );
        end
    end
end