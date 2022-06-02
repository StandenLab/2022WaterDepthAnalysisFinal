% function [p Rcrit] = moorelookup(n, R)
%
% This function looks up the 2-tailed p-value from the table in Zar (1999)
% Appendix B for the "Critical Values of R' for the Moore Test of Circular
% Uniformity" (Table B.39). Note that the alpha level has been doubled here
% to make the test two sided (as in the examples Zar presents for
% determining directionality of a non-normal distribution of rectangular
% coordinates).
%
% INPUT
% n: the length of your sample.
% R: R test statistic for the Moore Test of Circular Uniformity.
%
% OUTPUT
% p: p-value from the test.
% Rcrit: critical value of R.
%
% REQUIRES
% mooretable
%
% REFERENCES
% Zar, J.H. (1999). Biostatistical Analysis. Fourth edition. Upper Saddle
%   River, New Jersey: Prentice-Hall Inc.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [November 22, 2018]
%
% Edits:
% [20220428 | KL] - clean up function and comments
%__________________________________________________________________________

function [p, Rcrit] = moorelookup(n, R)

load('mooretable.mat','mooretable');
alpha = [0.999,0.2,0.1,0.05,0.02,0.01,0.002];
nn = mooretable(:,1);

% Find correct row of the table.
[easy, row] = ismember(n, nn);
if ~easy
   % Find closest value if no entry is present).
   row = length(nn) - sum(n<nn); 
   if row == 0
       error('N too small.');
   else
      warning('N=%d not found in table, using closest N=%d present.',n,nn(row)) %#ok<WNTAG>
   end
end

% Find minimal p-value and test-statistic.
idx = find(mooretable(row,2:end)<R,1,'last');
if ~isempty(idx)
  p = alpha(idx);
else
  p = 1;
end

% Lookup the critical value.
Rcrit = mooretable(row,idx+1);
if isempty(Rcrit)
    Rcrit=NaN;
end
end

