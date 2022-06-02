% function [mu2,r2]=circ_mean2(data,r)
%
% Calculates the second order mean based on Zar (1999) pg 600, 609-610.
% User inputs a set of mean angles and their associated r values. Data must
% be a single row.
%
% INPUT
% data: input angles in radians. Single row. Data across columns.
% r: input concentration parameters for angles. Single row. Data across
%   columns.
%
% OUTPUT
% mu2: second order angular mean.
% r2: concentration parameter for second order angular mean.
%
% REQUIRES
% circ_CorrAng
% mustHaveSameNaN
%
% REFERENCES
% Zar, J.H. (1999). Biostatistical Analysis. Fourth edition. Upper Saddle
%   River, New Jersey: Prentice-Hall Inc.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [November 21, 2018]
%
% Edits:
% [20220428 | KL] - cleaning function and comments
%__________________________________________________________________________

function [mu2,r2]=circ_mean2(data,r)

arguments
    data (1,:) double {mustBeVector}
    r (1,:) double {mustBeVector, mustHaveSameNaN(data,r)}
end

% Check for NaNs and remove.
tempdata=rmmissing(data);
tempr=rmmissing(r);

% Calculate rectangular coordinates for each mean angle.
X=nan(1,length(data));              % preset for speed
Y=nan(1,length(data));              % preset for speed
for i=1:size(tempdata,2)
    X(i)=tempr(i)*cos(tempdata(i));
    Y(i)=tempr(i)*sin(tempdata(i));
end

% Sum rectangular coordinates for each mean angle.
sumX=sum(X);
sumY=sum(Y);

% Calculate mean rectangular coordinates.
Xbar=sumX/size(tempdata,2);
Ybar=sumY/size(tempdata,2);

% Calculate r for set of mean angles.
r2=sqrt(Xbar^2+Ybar^2);

% Calculate the angle and place in the correct quadrant.
[mu2]=circ_CorrAng(Xbar/r2,Ybar/r2);

end