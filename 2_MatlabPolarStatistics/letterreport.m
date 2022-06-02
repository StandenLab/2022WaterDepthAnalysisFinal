% function [output]=letterreport(pvals,comps1,comps2,opts)
%
% Implementation of the algorithm written by Piepho (2004) for determining
% pairwise letter reports of p-values.
%
% INPUT
% pvals: a vector containing p-values
% comps1: a vector containing the group number for the first group in each
%   pariwise comparison.
% comps2: a vector containing the group number for the second group in each
%   pairwise comparison.
%
% NAME-VALUE ARGUMENTS
% pval: adjust the cutoff used for determing a comparison is statistically
%   significant. Default is 0.05.
%
% OUTPUT
% output: a string vector of letter comparisons. Ordered from group "1" of
%   the input to the last group of the input.
%
% REFERENCES
% Piepho, H.P. (2004). An algorithm for a letter-based representation of
% all-pairwise comparisons. Journal of Computational and Graphical
% Statistics, 13(2), 456-466.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [December 16, 2021]
%
% Edits:
% [20220503 | KL] - comments and cleaning up code
%__________________________________________________________________________

function [output]=letterreport(pvals,comps1,comps2,opts)
arguments
    pvals (1,:) double
    comps1 (1,:) double {mustBeEqualSize(pvals,comps1)}
    comps2 (1,:) double {mustBeEqualSize(pvals,comps2)}
    opts.pval (1,1) double = 0.05
end

%Pull all the groups.
groups=sort(unique([unique(comps1) unique(comps2)]));

%Initiate letterreport matrix & it's maximum value
matrix=ones(length(groups),1);
matmax=1;

%For each p-value
for i=1:length(pvals)
    if pvals(i)<opts.pval                                                   % if there is a statistically significant difference
        for m=1:size(matrix,2)                                              % iterate across columns of your putative letterreport matrix
            if ~isnan(matrix(comps1(i),m)) && ~isnan(matrix(comps2(i),m))   % if the current column connects the significant difference
                matmax=matmax+1;                                            % update maximum value
                matrix(:,end+1)=matrix(:,m);                                % duplicate the connecting column at the end of the matrix
                matrix(~isnan(matrix(:,end)),end)=matmax;                   % update duplicate column to have maxmax as non-nan values
                matrix(comps1(i),m)=NaN;                                    % get rid of first connecting row
                matrix(comps2(i),end)=NaN;                                  % get rid of second connecting row
                if m~=1                                                     % we keep the first row as is cause it seems necessary
                    [matrix]=absorb(matrix,m);                              % absorb the row if possible
                end
                [matrix]=absorb(matrix,size(matrix,2));                     % absorb the row if possible
            end
        end
    end
end
clearvars -except matrix pvals comps1 comps2 opts
matrix=matrix(:,sum(isnan(matrix))~=size(matrix,1)); % get rid of any all-nan columns

%Rearrange the columns so it's easier to assign letters in a logical
%fashion
rearr=nan(size(matrix));                    % preset rearr matrix
cnt=1;                                      % set count
for i=1:size(matrix,1)                      % for each row of your matrix
    %Determine if there's data in this row
    col=~isnan(matrix(i,:));
    findcol=find(col);
    sumcol=sum(col);    
    if sumcol == 1                          % if only one column has data
        rearr(:,cnt)=matrix(:,findcol);
        matrix(:,findcol)=NaN;
        cnt=cnt+1;
        if all(all(isnan(matrix)))
            break
        end
    elseif sumcol>1                         % if more than one column has data
        
        %Find the rows with non-nan values in matrix
        realrows=cell(1,length(findcol));
        for f=1:length(findcol)
            realrows{f}=find(~isnan(matrix(:,findcol(f))));
        end
        
        %Find the intersection of those sets
        for f=1:length(realrows)-1
            if f==1
                common=intersect(realrows{f},realrows{f+1});
            else
                common=intersect(common,realrows{f+1});
            end            
        end
        
        %Get rid of rows common to all sets
        for f=1:length(realrows)
            for c=1:length(common)
                realrows{f}(realrows{f}==common(c))=NaN;
            end
        end
        
        %Find the minimum row reference in each set
        minrow=cellfun(@min,realrows);
        
        %Order findcol as you like
        [~,idx]=sort(minrow);
        orderedcols=findcol(idx);
        
        %Add the necessary columns to rearr
        for j=1:length(orderedcols)
            rearr(:,cnt)=matrix(:,orderedcols(j));
            matrix(:,orderedcols(j))=NaN;
            cnt=cnt+1;
        end
        
        %Get out of the loop when matrix is only NaNs
        if all(all(isnan(matrix)))
            break
        end
    elseif sumcol==0 % if no columns have data
    else
        error('Something has gone wrong. The code requires further debugging. Contact Keegan Lutek (kklutek@gmail.com) or have fun debugging yourself.')
    end
end

%Re-number each column sequntially
for j=1:size(rearr,2)
    rearr(~isnan(rearr(:,j)),j)=j;
end
clearvars -except rearr pvals comps1 comps2 opts

%Turn that set of numbers into letters
cellmat=compose('%g',rearr);
count='a';
for j=1:size(cellmat,2)
    for i=1:size(cellmat,1)
        if ~isequal(cellmat{i,j},'NaN')
            cellmat{i,j}=count;
        else
            cellmat{i,j}='';
        end
    end
    count=char(count+1);
end

%Combine individual letters into strings for output.
output=strings(1,size(cellmat,1));
for i=1:size(cellmat,1)
    output(i)=string(cat(2,cellmat{i,:}));
end
end

% Custom validation function for the size of input arguments
function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Inputs must have equal size.';
    throwAsCaller(MException(eid,msg))
end
end

%Algorithm for the absorb step. "matrix" is the putative letterreport
%matrix, curr is the current working column.
function [matrix]=absorb(matrix,curr)
arguments
    matrix (:,:) double
    curr (1,1) double
end

adjcol1=find(~isnan(matrix(:,curr))); %find non-nan rows in the current column

%pull columns other than the current working column
checkcols=1:size(matrix,2);
checkcols(curr)=NaN;
checkcols=rmmissing(checkcols);

%iterate across the necessary columns
for p=checkcols
    thiscol=find(~isnan(matrix(:,p))); %find non-nan rows in pth column
    %if curr column is a subset of pth column, set curr column to NaNs
    if all(ismember(adjcol1,thiscol))
        matrix(:,curr)=NaN;
    end
end
end