% function [out]=VMdistRayleigh(data,r)
%
% This function checks that your data are normally distributed and
% directional. Designed to work for BOTH first order and second order sets
% of angles (or mean angles). If your dataset is first order, feed in only
% your angles. If your dataset is second order, feed in your mean angles
% and their associated 'r' values. All tests are based on Zar (1999) and
% references therein. See each function for specific references. INPUT and
% OUTPUT in square brackets are optional.
%
% INPUT
% data: a set of angles (or mean angles). Single row. Data across columns.
% [r]: a set of concentration coefficients if your data set is second
%   order. Single row. Data across columns.
%
% OUTPUT
% mu: mean angle
% r: mean reesultant vector length
% CI: 95% confidence interval for a first order set of angles.
% kappa: kappa concentration parameter for describing the von Mises'
%   distribution closest to your data.
% vMdP: p-value testing Ho: the distribution of your data matches a von
%   Mises' distribution.
% angvar: angular variance of your data
% R: Raleigh's test statistic
% Rp: p-value from Raleigh's test. Ho: the population is uniformly
%   distributed about the circle.
% [mu2]: second order mean angle.
% [r2]: resultant vector length for a set of mean angles.
% [LL]: 95% lower confidence limit for a set of mean angles. Note that
%   confidence limits about a second order mean angle are not necessarily
%   symmetrical.
% [UL]: 95% upper confidence limit for a set of mean angles. Note that
%   confidence limits about a second order mean angle are not necessarily
%   symmetrical.
% [Rp]: p-value resulting from a parametric or non-parametric test for
%   directionality in a set of mean angles. If data are decided to be normal,
%   the test is run according to Hotelling (1931) and an F test statistic is
%   calculated. Ho: There is no mean population direction. If the data are
%   non-normal, a Moore-Raleigh test is run and an R' test statistic is
%   calculated. Ho: The population from which the sample of means came is
%   uniformly distributed around the circle.
% [R]: Test statistic for either of the tests mentioned in [Rp].
% [Rcrit]: Critical value for either of the test mentioned in [Rp].
%
% REQUIRES
% AngStats
% circ_kappa (Berens, 2009)
% circ_kuipertest (Berens, 2009)
% circ_CL2
% circ_mean (Berens, 2009)
% circ_mean2
% circ_rayleigh2
% circ_vmrnd (Berens, 2009)
%
% REFERENCES
% Berens, P. (2009). CircStat: A Matlab Toolbox for Circular Statistics.
%   Journal of Statistical Software. 31(10). [code available on Matlab file
%   exchange]
% Zar, J.H. (1999). Biostatistical Analysis. Fourth edition. Upper Saddle
%   River, New Jersey: Prentice-Hall Inc.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [October 1, 2018]
%
% Edits:
% [20220426 | KL] - cleaning function and comments
%__________________________________________________________________________

function [out]=VMdistRayleigh(data,r)

arguments
    data (1,:) double
    r (1,:) double {mustHaveSameNaN(data,r)} = data
end

switch nargin
    case 1
        % Get rid of any NaNs in the dataset.
        tempdata=rmmissing(data);
        
        if ~isempty(tempdata)
            [mu]=circ_mean(tempdata');                      % calculate circular mean of your input data
            [kappa]=circ_kappa(tempdata');                  % calculate kappa parameter of the von Mises distribution for your input data
            [vonmis]=circ_vmrnd(mu,kappa,length(tempdata)); % simulate a random sample the same length as your input data with parameters mu and kappa
            
            % If your data is above the required sample size, check whether
            % your data follows a von Mises distribution
            kuiper=nan(1,3);
            if size(tempdata,2) > 4
                [kuiper(1),kuiper(2), kuiper(3)]=kuiperTest(tempdata',vonmis);
            else
                warning("Sample size not large enough for Kuiper Test. Must be > 4.")
            end
            
            % If your data is approximately von Mises distributed, run a
            % Rayleigh test, calculate confidence interval and angular
            % variance.
            if kuiper(1) > 0.049
                [raleigh]=AngStats(tempdata');
            end
            
            % Assign for output
            out.mu=mu;                      % angular mean
            out.kappa=kappa;                % kappa parameter
            out.vMdP=kuiper(1);             % p-value for von Mises distribution test | H0: your data is approximately von Mises distributed
            if exist('raleigh','var')       % if your data has a von Mises distribution, assign output
                out.CI = raleigh.CI;        % confidence interval
                out.r=raleigh.r;            % r - "concentration parameter"
                out.angvar=raleigh.angvar;  % angular variance
                out.R=raleigh.R;            % Rayleigh's test statistic
                out.Rp=raleigh.Rp;          % Rayleigh's p-value
            else                            % if your data isn't von Mises distributed, assign NaNs
                out.CI=NaN;
                out.r=NaN;
                out.angvar=NaN;
                out.R=NaN;
                out.Rp=NaN;
                
            end
        else %Data is all NaNs
            warning("Input data is empty.")
            out.mu=NaN;
            out.kappa=NaN;
            out.vMdP=NaN;
            out.CI=NaN;
            out.r=NaN;
            out.angvar=NaN;
            out.R=NaN;
            out.Rp=NaN;
        end
    case 2 %Two input arguments
        
        % Get rid of any NaNs in the dataset.
        tempdata=rmmissing(data);
        tempr=rmmissing(r);
        
        % Calculate second order angular mean.
        [mu,r]=circ_mean2(tempdata,tempr);
        
        % If data length is sufficient, calculate second order confidence
        % interval and run second order Rayleigh's test (parametric or
        % non-parametric).
        if size(tempdata,2) > 2
            [Ralp,Ralstat,Ralcrit,normal] = circ_rayleigh2(tempdata,tempr);     % second order Rayleigh's test
            if Ralp < 0.05                                                      % if the data is directional, calculate confidence limits
                [LL,UL]=circ_CL2(tempdata,tempr,mu);                            % second order angular mean
            else                                                                % the data is not directional, cannot calculate confidence limits
                LL='data not directional. cannot compute confidence limits.';
                UL='data not directional. cannot compute confidence limits.';
            end
        else % insufficient data to run statistical tests.
            LL='n too small';
            UL='n too small';
            Ralp='n too small';
            Ralstat='n too small';
            Ralcrit='n too small';
            normal='n too small';
        end
        
        %Assign for output
        out.mu2=mu;
        out.r2=r;
        out.LL=LL;
        out.UL=UL;
        out.Rp=Ralp;
        out.R=Ralstat;
        out.Rcrit=Ralcrit;
        out.normal=normal;
end
end