function output = callsdpagmp_neos(interfacedata)

% CALLSDPAGMP_NEOS.m Send SDPA-GMP problem to NEOS solver from YALMIP 

% ----------------------------------------------------------------------- %
%
%     Copyright (C) 2017 Hugo Tadashi
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
% ----------------------------------------------------------------------- %

% Check for xmlrpc-client-1.1.1.jar
if ~exist('xmlrpc-client-1.1.1.jar','file')
    error('Please add xmlrpc-client-1.1.1.jar to your dynamic class path');
end

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Check if param.sdpa file exists in pwd
cleanup = 0;
if ~exist([pwd,filesep,'param.sdpa'],'file')
    cleanup = 1;
    % write it with default parameters, otherwise failure!
    fID = fopen([pwd,filesep,'param.sdpa'],'w');
    fprintf(fID,'200         unsigned int    maxIteration;                 \n');
    fprintf(fID,'1.0E-25     double          0.0 < epsilonStar;            \n');
    fprintf(fID,'1.0E6       double          0.0 < lambdaStar;             \n');
    fprintf(fID,'2.0         double          1.0 < omegaStar;              \n');
    fprintf(fID,'-1.0E25     double          lowerBound;                   \n');
    fprintf(fID,'1.0E25      double          upperBound;                   \n');
    fprintf(fID,'0.1         double          0.0 <= betaStar <  1.0;       \n');
    fprintf(fID,'0.2         double          0.0 <= betaBar  <  1.0, betaStar <= betaBar;\n');
    fprintf(fID,'0.7         double          0.0 < gammaStar <  1.0;       \n');
    fprintf(fID,'1.0E-25     double          0.0 < epsilonDash;            \n');
    fprintf(fID,'200         precision;                                    \n');
    fclose(fID);
end

% Bounded variables converted to constraints
if ~isempty(ub)
    % addbounds was renamed somewhere along the way in YALMIP?
    try
        [F_struc,K] = addbounds(F_struc,K,ub,lb);
    catch
        [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
    end
end

% Convert from internal (sedumi) format
[mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);

if options.verbose==0
    options.sdpa.print = 'no';
else
    options.sdpa.print = 'display';
end

if options.savedebug
    ops = options.sdpa;
    save sdpadebug mDIM nBLOCK bLOCKsTRUCT c F ops
end

if options.showprogress
    showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL SDPA-GMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export to SDPA-GMP, solve and inport results
FF = cellfun(@full,F,'UniformOutput',false);
[nmax,mmax] = size(FF);
BS = abs(bLOCKsTRUCT);
for nn=1:nmax               % loop to construc SDPA problem
    for mm=1:mmax
        if isempty(FF{nn,mm})
            FF{nn,mm} = zeros(BS(nn));
        end
    end
end

tic
inputSDPA  = 'sdpagmp_in.dat-s';
outputSDPA = 'sdpagmp_out.out';
header     = 'Input from YALMIP';

gensdpagmpfile(inputSDPA,mDIM,nBLOCK,bLOCKsTRUCT,c,FF,header);          % write SDPA-GMP input file

% Send job to NEOS server and write solution to outputSDPA file
disp('Creating NEOS server interface')
disp('')
neos = NeosSDPAInterface();

queue = neos.get_queue();
disp('Job queue:')
disp('')
disp(queue)

precision = 'var';
xml_string = build_xml_string(inputSDPA,'param.sdpa',precision);

outputSDPA_filepath = [pwd filesep outputSDPA];
neos.submit_job(xml_string,outputSDPA_filepath);

% import result
[objVal,x,X,Y,INFO] = sdpagmp_read_output(outputSDPA,full(mDIM),full(nBLOCK),full(bLOCKsTRUCT));
solvertime = toc;

% Clean up tmp files created in this directory
delete(inputSDPA);
delete(outputSDPA);
if cleanup
    delete('param.sdpa');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here onwards, like in YALMIP native callsdpa
% Create variables in YALMIP internal format
Primal = x;

Dual = [];
for i = 1:length(Y)
    Dual = [Dual;Y{i}(:)];
end

Slack = [];
if options.saveduals
    for i = 1:length(X)
        Slack = [Slack;X{i}(:)];
    end
end

switch (INFO.phasevalue)
    case 'pdOPT'
        problem = 0;
    case {'noINFO','pFEAS','dFEAS'}
        problem = 3;
    case {'pdFEAS'}
        problem = 4;
    case 'pFEAS_dINF'
        problem = 2;
    case 'pINF_dFEAS'
        problem = 1;
    case 'pUNBD'
        problem = 2;
    case 'dUNBD'
        problem = 1;
    case 'pdINF'
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput
    solveroutput.objVal = objVal;
    solveroutput.x = x;
    solveroutput.X = X;
    solveroutput.Y = Y;
    solveroutput.INFO = INFO;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.mDIM = mDIM;
    solverinput.nBLOCK=nBLOCK;
    solverinput.bLOCKsTRUCT=bLOCKsTRUCT;
    solverinput.c=c;
    solverinput.F=F;
else
    solverinput = [];
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);

end

function xml_string = build_xml_string(sdpa_file,param_sdpa,precision)
% XML header
xml_string = sprintf('<document>\n<category>%s</category>\n<solver>%s</solver>\n<inputType>%s</inputType>\n','sdp','SDPA','SPARSE_SDPA');

% File in sparse SDPA format
dat = fileread(sdpa_file);
insert = sprintf('<dat><![CDATA[%s]]></dat>\n',dat);
xml_string = strcat(xml_string,insert);

% MATLAB binary file (not used)
insert = sprintf('<mat><![CDATA[]]></mat>\n');
xml_string = strcat(xml_string,insert);

% param.sdpa file
param = fileread(param_sdpa);
insert = sprintf('<param><![CDATA[%s]]></param>\n',param);
xml_string = strcat(xml_string,insert);

% Precision
insert = sprintf('<PRECISION><![CDATA[%s]]></PRECISION>\n',precision);
xml_string = strcat(xml_string,insert);

% Comments (not used)
insert = sprintf('<comment><![CDATA[]]></comment>\n');
xml_string = strcat(xml_string,insert);

% End of xml
xml_string = strcat(xml_string,'</document>');

% Parse XML
xml_string = xmlwrite(xmlreadstring(xml_string));

end















