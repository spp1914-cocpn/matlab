classdef SimulationInfoPrinter < handle
    % Simple class to print the progress of a running (Monte Carlo)
    % simulation on the screen (i.e., the command window):
    % Out format is: simulationName: numRun/maxRuns | numTimeStep/timeStepsPerRun (percentage)
    % Example output: 'My simulation: 5/10 | 500/1000 (45%)'
    
    properties (GetAccess = private, SetAccess = immutable)
        numRuns;
        numTimeStepsPerRun;
        printFormatString;
    end
    
    properties (Access = private)
        printingEnabled;
        lastPrintedLength = 0;
    end
    
    properties (Dependent, Access = private)
        total;
    end
    
    methods
        function total = get.total(this)
            total = this.numRuns * this.numTimeStepsPerRun;
        end
    end
    
    methods (Access = public)
        function this = SimulationInfoPrinter(simulationName, numRuns, numTimeStepsPerRun)
            this.numRuns = numRuns;
            powersOfTen = 10 .^ [1:100];
            idx = find(numRuns ./ powersOfTen, 1);
            this.numTimeStepsPerRun = numTimeStepsPerRun;
            idx2 = find(numTimeStepsPerRun ./ powersOfTen, 1);
            this.printFormatString = sprintf('%s: %%0%dd/%%%dd | %%0%dd/%%%dd (%%3.0f%%%%)', ...
                   simulationName, idx, idx, idx2, idx2);
               
           this.turnOn();
        end
        
        function turnOn(this)
            this.printingEnabled = true;
        end
        
        function turnOff(this)
            this.printingEnabled = false;
        end
        
        function printSimulationStart(this)
            if this.printingEnabled
                this.lastPrintedLength = fprintf(this.printFormatString, 0, this.numRuns, ...
                    0, this.numTimeStepsPerRun, 0);
            end
        end
        
        function printProgress(this, numRun, numTimeStep)
            currentRun = max(numRun - 1, 0) * this.numTimeStepsPerRun + numTimeStep;
            if this.printingEnabled
                progress = floor(100 * currentRun / this.total); % percentage
                delims = repmat('\b', 1, this.lastPrintedLength); %\b requires 2 byte each
                this.lastPrintedLength = fprintf([delims this.printFormatString], ...
                    numRun, this.numRuns, numTimeStep, this.numTimeStepsPerRun, progress) - numel(delims) / 2;
            end
        end
        
        function printSimulationEnd(this)
            if this.printingEnabled
                fprintf('\n');
            end
        end
    end
end

