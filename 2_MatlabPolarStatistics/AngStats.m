% function [out]=AngStats(data)
%
% Runs Rayleigh's test for a unimodal distribution and calculates the
% angular mean, angular variance and confidence interval for a first order
% circular dataset.
%
% INPUT
% data: input data. Data across rows. One column.
%
% OUTPUT
% data: your input data.
% mu: angular mean.
% k: length of your input data.
% r: concentration parameter
% angvar: angular variance.
% CI: confidence interval.
% R: Rayleigh's test statistic. H0 - data is uniformly distributed.
% Rp: Rayleigh's p-value.
% output: the above output in a table.
%
% REQUIRES
% circ_mean (Berens, 2009)
%
% REFERENCES
% Berens, P. (2009). CircStat: A Matlab Toolbox for Circular Statistics.
%   Journal of Statistical Software. 31(10). [code available on Matlab file
%   exchange]
%
%__________________________________________________________________________
% Written by: Emily M. Standen [October 1, 2007]
%
% Edits:
% [20220427 | KL] - updates to comments and cleaning function.
%__________________________________________________________________________

function [out]=AngStats(data)

arguments
    data (:,1) double {mustBeNumeric}
end

% Remove any NaNs from the data
datashort=rmmissing(data);

if isempty(datashort)   % if input data is empty
    
    % Assign output.
    out.data=data;
    out.mu=NaN;
    out.k=NaN;
    out.r=NaN;
    out.angvar=NaN;
    out.CI=NaN;
    out.R=NaN;
    out.Rp=NaN;
    out.output=NaN;
else                    % if there are finite values in input data
    
    % Get the angular mean
    mu = circ_mean(datashort);
    
    % Get the dispersion of the angles
    y=nan(size(datashort));
    x=nan(size(datashort));
    for i=1:size(datashort,1)
        y(i)=sin(datashort(i));
        x(i)=cos(datashort(i));
    end
    k=sum(isfinite(y));
    Y=(sum(y))/k;
    X=(sum(x))/k;
    
    % Calculate concentration parameter
    r=sqrt(Y^2 +X^2);
    
    % Calculate angular variance
    angvar=2*(1-r); % Zar pg.604 in units of radians squared
    
    % Calculate Rayleigh's test statistic and confidence intervals
    R = r*k;
    if r<=0.9 && r>sqrt(5.024/(2*k))
        CI = acos(sqrt(2*k*(2*(R)^2 - k*5.024)/(4*k - 5.024))/(R));
    elseif r>= 0.9
        CI = acos(sqrt(k^2 - (k^2 - (R)^2)*(2.71828^(5.024/k)))/(R));
    else
        CI="r is too small";
    end
    
    % Calculate Rayleigh's p-value.
    Rp = exp(sqrt(1 + 4*k + 4*(k^2 - R^2)) - (1+2*k));
    
    % Generate output tables.
    if validateCI(CI)   %check is CI is the above string
        output=table(mu,k,NaN,angvar,NaN,R,Rp);                                     % radian output table
        output.Properties.VariableNames=["mu","k","r","angvar","CI","R","Rp"];      % update variable names
        degoutput=table(rad2deg(mu),NaN,NaN,NaN,NaN,NaN,NaN);                       % degrees output table
        degoutput.Properties.VariableNames=["mu","k","r","angvar","CI","R","Rp"];   % update variable names
        output=[output;degoutput];                                                  % concatenate together
        output.Properties.RowNames=["radians","degrees"];                           % add row names
    else
        output=table(mu,k,r,angvar,CI,R,Rp);                                        % radian output table
        output.Properties.VariableNames=["mu","k","r","angvar","CI","R","Rp"];      % update variable names
        degoutput=table(rad2deg(mu),NaN,NaN,rad2deg(angvar),NaN,NaN,NaN);           % degrees output table
        degoutput.Properties.VariableNames=["mu","k","r","angvar","CI","R","Rp"];   % update variable names
        output=[output;degoutput];                                                  % concatenate together
        output.Properties.RowNames=["radians","degrees"];                           % add row names
    end
    
    % Assign output
    out.data=data;
    out.mu=mu;
    out.k=k;
    out.r=r;
    out.angvar=angvar;
    out.CI=CI;
    out.R=R;
    out.Rp=Rp;
    out.output=output;
end
end

% Validation function for the variable CI.
function out=validateCI(in)
    out=0;
    if isstring(in)
        if matches(in,"r is too small")
            out=1;
        end
    end
end