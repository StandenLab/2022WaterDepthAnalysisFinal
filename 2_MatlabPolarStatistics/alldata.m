% function [dat,treat,dif]=alldata(structin,varname,opts)
%
% Collects and concatenates data from a nested structure with data split by
% treatments. e.g. structin.TREAT.varname.
%
% INPUT
% structin: input structure. First level of nesting should be a set of
%   structures with data for each condition of interest.
% varname: string with your variable name of interest.
%
% NAME-VALUE INPUTS
% angdif: compute the angular difference or not. 0 = don't compute. 1 = do
%   compute. Default is 0.
%
% OUTPUT
% dat: data concatenated into a single cell vector. Each cell is a
%   position on the fish.
% treat: cell vector with each cell containing a set of strings. Each
%   individual string denotes the treatment for each value in dat.
% dif: if asked for, the angular difference. Organized as for dat.
%
% REQUIRES
% circ_dist (Berens, 2009)
% circ_mean (Berens, 2009)
%
% REFERENCES
% Berens, P. (2009). CircStat: A Matlab Toolbox for Circular Statistics.
% Journal of Statistical Software. 31(10). [code available on Matlab file
% exchange]
%
%__________________________________________________________________________
% Written by: Keegan Lutek [December 23, 2021]
%__________________________________________________________________________

function [dat,treat,dif]=alldata(structin,varname,opts)
arguments
    structin (:,:) struct
    varname (:,:) string
    opts.angdif (:,:) = 0
end
conds=fieldnames(structin);
cnt=0;
for c=1:length(conds)                                               % for all conditions
    if isfield(structin.(conds{c}),varname)                         % if the variable exists for the current condition
        cnt=cnt+1;
        working=structin.(conds{c});
        
        % Convert to cell if necessary
        if ~iscell(working.(varname))
            working.(varname)=num2cell(working.(varname),[1 2]);
            warning('%s converted to single cell',varname);
        end
        
        % Initialize variables when necessary
        if cnt==1
            dat=cell(1,length(working.(varname)));
            treat=cell(1,length(working.(varname)));
            if opts.angdif==1
                dif=cell(1,length(working.(varname)));
            end
        end
        % Add data to the end of variables
        for j=1:length(working.(varname))
            dat{j}=[dat{j} working.(varname){j}];
            treat{j}=[treat{j} repelem(c,(length(working.(varname){j})))];
            if opts.angdif==1
                thesedif=abs(circ_dist(working.(varname){j},circ_mean(working.(varname){j})));
                dif{j}=[dif{j} thesedif];
                clear thesedif
            end
        end
    end
end
if opts.angdif==0
    dif=NaN;
end
end
