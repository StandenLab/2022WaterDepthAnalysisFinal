% function [p]=HermansRasson2P(sample)
%
% Runs the Hermans-Rasson test as defined by Landler et al. (2019) to test
% for a non-uniform circular distribution. This can be used as an
% alternative to Rayleigh's test for data that does not have a unimodal von
% Mises distribution.
%
% H0: data is randomly distributed about the circle.
%
% INPUT
% sample: your sample data in RADIANS.
%
% OUTPUT
% p: the p-value for your Hermans-Rasson test.
%
% REFERENCES
% Landler, L., Ruxton, G.D., and Malkemper, E.P. (2019). The Hermans-Rasson
% test as a powerful alternative to the Rayleigh test for circular
% statistics in biology. 19, 30.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [July 12, 2021]
%
% Edits:
% [20220426 | KL] - clean up comments and documentation
%__________________________________________________________________________

function [p]=HermansRasson2P(sample)

arguments
    sample (1,:) {mustBeNumeric}
end

    % Generate required presets
    univals=9999;           % number of iterations for calculating p-value
    n=length(sample);       % find the length of your sample
    testset=nan(1,univals); % preset for speed.
    
    % Run iterations to estimate the test statistic for a randomly
    % distributed sample of the same length as your input data.
    for f=1:univals
        theta = rand(1,n)*(2*pi);           %Generate a randomly distributed sample about a circle of the same length as your input data
        testset(f)=HermansRasson2T(theta);  %Test statistic for the random sample.
    end
    
    %Find the test statistic for your sample
    Tsample=HermansRasson2T(sample);
    
    % Calculate the p-value: 
    % Let "Q" be the number of random samples that give a test statistic
    % equal to or greater in magnitude to that of the original sample. Let
    % "m" be the number of random samples you have generated. The p-value
    % of the test is given by (Q + 1)/(m + 1). For details see Landler et
    % al. (2019).
    counter=1;
    for j=1:univals
        if testset(j)>=Tsample
            counter=counter+1;
        end
    end
    p=counter/(univals+1); % Calculate your p-value.
end

% Calculates the new formulation of the Hermans-Rasson test statistic as
% defined by Landler et al. 2019 BMC Ecology
function [T]=HermansRasson2T(sample)
arguments
    sample (1,:) {mustBeNumeric}
end
    n=length(sample); % find the length of your sample
    
    % Calculate the numerator of the test statistic
    total=0;    % initialize the sum
    for i=1:n
        for j=1:n
            total=total + abs(abs(sample(i)-sample(j))-pi)-(pi/2);
            total=total - (2.895*(abs(sin(sample(i)-sample(j)))-(2/pi)));
        end
    end
    
    % Calculate the test statistic
    T=total/n;
end