% function [output]=pairwiseletters(tablein)
%
% Essentially a wrapper for letterreport that allows the user to iterate
% through a series of body positions.
%
% INPUT
% tablein: input data in a table. Row names must be each comparison in the
%   following form "1 - 2".
%
% OUTPUT
% output: a string array of connecting letters. Each string represents a
%   condition. Order should be the order they appear in tablein from top to
%   bottom. Each column in your cell array is a column in your input table.
%
% REQUIRES
% letterreport
%
%__________________________________________________________________________
% Written by: Keegan Lutek [December 23, 2021]
%__________________________________________________________________________

function [output]=pairwiseletters(tablein)
arguments
    tablein (:,:) table
end

% Determine the first and second conditions for each comparison based on
% the rownames of the input table.
for i=1:size(tablein,1)
    splitcomp=strsplit(tablein.Properties.RowNames{i},' - ');
    comps1(i)=string(splitcomp{1});
    comps2(i)=string(splitcomp{2});
end

% Initiate output table
output=cell(1,size(tablein,2));

% Generate a letter report for the pairwise comparisons
for j=1:size(tablein,2)
    pvals=table2array(tablein(:,j));
    if all(~isnan(pvals))
    [output{j}]=letterreport(pvals,comps1,comps2);
    end
end
end