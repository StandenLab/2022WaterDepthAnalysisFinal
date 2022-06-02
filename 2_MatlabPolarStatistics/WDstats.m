% function [outtable,data]=WDstats(datain,cols)
%
% Essentially a wrapper function for HermansRasson2P and VMdistRayleigh that
% runs the inital stats for the WaterDepth manuscript polar variables. Data
% are first tested for a non-uniform direction using the HR test (after
% Landler et al. 2019). If the data is non-uniformly distributed, the data
% is tested for a unimodal distribution using Rayleigh's test.
%
% INPUT
% datain: a cell vector of all of your data points. Each cell is a
%   condition.
% cols: a string vector that contains labels for each of the cells in
%   datain.
%
% OUTPUT
% outtable: results of the above tests in table form.
%   "HRp" is the Hermans-Rasson p-value [p<0.05 = non-uniform distribution]
%   "VMp" is the von Mises p-value [p>0.05 = von Mises distributed data]
%   "Rp" is the Rayleigh's test p-value [p<0.05 = unimodal, directional]
%   "mu" is the angular mean of your data
%   "r" is the mean resultant vector length
%   "angvar" is the angular variance
%   "CI" is the 95% confidence interval for each mean
% data: datain returned to the user.
%
% REQUIRES
% AngStats
% circ_CL2
% circ_kappa (Berens, 2009)
% circ_kuipertest (Berens, 2009)
% circ_mean (Berens, 2009)
% circ_mean2
% circ_raleigh2
% circ_vmrnd (Berens, 2009)
% HermansRasson2P
% VMdistRayleigh
%
% REFERENCES
% Berens, P. (2009). CircStat: A Matlab Toolbox for Circular Statistics.
% Journal of Statistical Software. 31(10). [code available on Matlab file
% exchange]
%
%__________________________________________________________________________
% Written by: Keegan Lutek [December 23, 2021]
%__________________________________________________________________________

function [outtable,data]=WDstats(datain,cols)
arguments
    datain (1,:) cell
    cols (1,:) string
end

% Preset output size
HRp=nan(1,size(datain,2));
mu=nan(size(HRp));
r=nan(size(HRp));
angvar=nan(size(HRp));
CI=nan(size(HRp));
Rp=nan(size(HRp));
VMp=nan(size(HRp));

for j=1:size(datain,2) % For each cell in datain
    % Run the Hermans-Rasson test
    HRp(1,j)=HermansRasson2P(datain{j});
    % If data is not uniformly distributed about the circle, run Rayleigh's
    if HRp(j) < 0.05
        Rayleigh=VMdistRayleigh(datain{j});
        % If data is unimodal, calculate mean etc.
        if Rayleigh.Rp < 0.05
            mu(1,j)=Rayleigh.mu;
            r(1,j)=Rayleigh.r;
            angvar(1,j)=Rayleigh.angvar;
            CI(1,j)=Rayleigh.CI;
        else % Not unimodal, cannot calculate mean etc.
            mu(1,j)=NaN;
            r(1,j)=NaN;
            angvar(1,j)=NaN;
            CI(1,j)=NaN;
        end
        Rp(1,j)=Rayleigh.Rp;
        VMp(1,j)=Rayleigh.vMdP;
    else % Uniformly distributed about cicle, cannot run Rayleighs etc.
        Rp(1,j)=NaN;
        VMp(1,j)=NaN;
        mu(1,j)=NaN;
        r(1,j)=NaN;
        angvar(1,j)=NaN;
        CI(1,j)=NaN;
    end
    
    clear Rayleigh
end

% Output table
outtable=splitvars(table([HRp;VMp;Rp;mu;r;angvar;CI]'));
outtable.Properties.RowNames=cols;
outtable.Properties.VariableNames=["HRp","VMp","Rp","mu","r","angvar","CI"];

% Output original data
data=cell(1,length(datain));
for j=1:length(datain)
    data{j}=rmmissing(datain{j});
end
end % end of main function
