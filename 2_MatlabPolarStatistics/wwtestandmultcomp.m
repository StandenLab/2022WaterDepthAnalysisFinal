% function [outF,outP,outP_bonf]=wwtestandmultcomp(datain,condin,positions)
%
% Essentially a wrapper for circ_wwtest. Allows the user to iterate across
% positions on the fish and run pairwise comparisons if required.
%
% INPUT
% datain: input data in a cell vector. Each cell is a body position.
% condin: input conditions in a cell vector. Each cell is a body position.
%   The length of the double in each cell must match those of datain.
%   Conditions should be numbers.
% positions: string input denoting the positions in each cell of datain.
%
% OUTPUT
% outF: a table with the results of your Watson-Williams test.
% outP: pairwise comparisons of the mean of each group if required following
%   the Watson-Williams test; p-values are uncorrected.
% outP_bonf: pairwise comparisons of the mean of each group if required
%   following the Watson-Williams test; p-values are bonferroni corrected.
%
% REQUIRES
% circ_wwtest (Berens, 2009)
%
% REFERENCES
% Berens, P. (2009). CircStat: A Matlab Toolbox for Circular Statistics.
%   Journal of Statistical Software. 31(10). [code available on Matlab file
%   exchange]
%
%__________________________________________________________________________
% Written by: Keegan Lutek [December 23, 2021]
%__________________________________________________________________________

function [outF,outP,outP_bonf]=wwtestandmultcomp(datain,condin,positions)
arguments
    datain (1,:) cell
    condin (1,:) cell
    positions (1,:) string
end

ncond=length(unique(condin{1})); %Find the number of conditions

%Preset matrix sizes
Fp=cell(1,length(datain));
FTable=cell(1,length(datain));
statmatrix=nan(length(datain),4);
multp=nan(ncond*(ncond-1)/2,1,length(datain));
multpcomp=cell(1,length(datain));

for j=1:length(datain)                                                          % for each position
    s = warning('error', 'MATLAB:oops');                                        % set the warning in circ_wwtest to an error so you can catch it
    try
    [Fp{j},FTable{j}]=circ_wwtest(datain{j},condin{j});                         % circ_wwtest across ALL groups
    catch                                                                       % if the above warning is thrown
        Fp{j}=NaN;                                                              % set p-value to NaN (test not applicable)
        FTable{j}=cell(4,6);                                                    % set Table data to NaN (test not applicable)
    end
    
    %Pull out the needed data from above
    if sum(sum(cellfun(@nanmean,FTable{j}(2:4,2:6)),'omitnan'),'omitnan') > 0   % check for whether there is data in the table (that's saved as a cell array) or not.
        statmatrix(j,1)=FTable{j}{2,5};
        statmatrix(j,2)=FTable{j}{2,2};
        statmatrix(j,3)=FTable{j}{4,2};
        statmatrix(j,4)=Fp{j};
    end
    
    % If there is a difference somewhere
    if Fp{j} < 0.05 
        cnt=1;                                                                                  % initialize count variable
        for i=1:ncond-1                                                                         % for (almost) all conditions
            for k=i+1:ncond                                                                     % generate the second part of the paired comparison
                theserefs=find(condin{j}==i | condin{j}==k);                                    % row references for relevant data
                try
                    [multp(cnt,j),~]=circ_wwtest(datain{j}(theserefs),condin{j}(theserefs));    % run two-sample wwtest
                catch                                                                           % catch error if needed (test not applicable)
                    multp(cnt,j)=NaN;
                end
                multpcomp{j}(cnt,1)=sprintf("%i - %i",i,k);                                     % generate each paired comparison marker
                cnt=cnt+1;
            end            
        end
        
    else
        multp(:,j)=nan(ncond*(ncond-1)/2,1);
        cnt=1;
        for i=1:ncond-1                                     %for (almost) all conditions
            for k=i+1:ncond                                 %generate the second part of the paired comparison
                multpcomp{j}(cnt,1)=sprintf("%i - %i",i,k); %generate each paired comparison marker
                cnt=cnt+1;
            end
        end
    end
    warning(s)
end

%Generate a table with the F test info from the across ALL groups
%comparison
outF=splitvars(table(statmatrix));
outF.Properties.RowNames=positions;
outF.Properties.VariableNames=["F","numdf","dendf","pval"];

%Generate a table with pairwise comparisons (uncorrected p values) as
%required
outP=splitvars(table(multp));
outP.Properties.RowNames=multpcomp{1};
outP.Properties.VariableNames=positions;

%Generate a table with (bonferroni corrected) pairwise comparisons as
%required
adjp=multp.*size(multp,1);
adjp(adjp>1)=1;
outP_bonf=splitvars(table(adjp));
outP_bonf.Properties.RowNames=multpcomp{1};
outP_bonf.Properties.VariableNames=positions;
end
