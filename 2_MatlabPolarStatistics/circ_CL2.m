% function [LL,UL]=circ_CL2(data,r,mu2)
%
% This function calculates the 95% confidence limits for a set of mean
% angles. Based on Zar (1999) pg 611-614. User inputs the set of mean
% angles, the concentration parameters for that set of angles and the
% second order mean of that set of angles. All inputs must be in a single
% row and in RADIANS.
%
% INPUT
% data: input angles in radians. Single row. Data across columns.
% r: input concentration parameters for angles. Single row. Data across
%   columns.
% mu2: second order angular mean for your dataset.
%
% OUPUT
% LL: lower limit of a 95% confidence interval. May be asymmetric about the
%   mean.
% UL: upper limit of a 95% confidence interval. May be asymmetric about the
%   mean.
%
% REQUIRES
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

function [LL,UL]=circ_CL2(data,r,mu2)

arguments
    data (1,:) double {mustBeVector}
    r (1,:) double {mustBeVector, mustHaveSameNaN(data,r)}
    mu2 (1,1) double
end

%% Initial data [k,X,Y,XY]
% Set # of samples (k)
k=size(data,2);

% Calculate rectangular coordinates for each mean angle.
X=nan(1,length(data));      % preset for speed
Y=nan(1,length(data));      % preset for speed
for i=1:length(data)
    X(i)=r(i)*cos(data(i));
    Y(i)=r(i)*sin(data(i));
end

% Calculate rectangular coodinates ^2 for each mean angle.
X2=X.^2;
Y2=Y.^2;

% Calculate product of rectangular coordinates for each mean angle.
XY=X.*Y;
%% Sums of things & mean rectangular coordinates
% Sum rectangular coordinates.
sumX=sum(X);
sumY=sum(Y);

% Sum rectangular coordinates ^2.
sumX2=sum(X2);
sumY2=sum(Y2);

% Sum product of rectangular coordinates.
sumXY=sum(XY);

% Calculate mean rectangular coordinates.
Xbar=sumX/k;
Ybar=sumY/k;
%% Sums of squares for X, Y and XY
% Calculate sum of squares.
ssX=sumX2-((sumX^2)/k);
ssY=sumY2-((sumY^2)/k);
ssXY=sumXY-((sumX*sumY)/k);
%% Calculate "A"
A=(k-1)/ssX;
%% Calculate "B"
B=-(((k-1)*ssXY)/(ssX*ssY));
%% Calculate "C"
C=(k-1)/ssY;
%% Find F test statistic
F=finv(0.95,2,(k-2));
%% Calculate "D"
D=(2*(k-1)*(1-((ssXY^2)/(ssX*ssY)))*F)/(k*(k-2));
%% Calculate "H"
H=A*C-B^2;
%% Calculate "G"
G=(A*(Xbar^2))+(2*B*Xbar*Ybar)+(C*(Ybar^2))-D;
%% Calculate "U"
U=(H*(Xbar^2))-(C*D);
%% Calculate "V"
V=sqrt(D*G*H);
%% Calculate "W"
W=(H*Xbar*Ybar)+(B*D);
%% Calculate "b1" and "b2" [upper and lower quantiles]
b1=(W+V)/U;
b2=(W-V)/U;
%% Calculate "M" & find the upper and lower limits.
% For "b1" [Upper Limit]
M1=sqrt(1+(b1^2));
% Note that the reference angle for cos(UL) = 1/M1 & sin(UL) = b1/M1

% Find the "correct" quadrant for angle denoted by acos(1/M1) & asin(b1/M1)
[UL]=circ_CorrAng(1/M1,b1/M1);

% Find the second option for the upper limits - see Zar pg 613
if UL+pi > 2*pi
    UL2=(UL+pi)-2*pi;
else
    UL2=UL+pi;
end

% Find the distances from UL/UL2 to mu2
% Forward UL:
if mu2 < UL
    fUL=abs(UL-mu2);
elseif mu2 > UL
    fUL=((2*pi)-mu2)+UL;
end
% Forward UL2:
if mu2 < UL2
    fUL2=abs(UL2-mu2);
elseif mu2 > UL2
    fUL2=((2*pi)-mu2)+UL2;
end

% Sort those distances to find the smallest 
tmp=[fUL,fUL2];
[~,rnkUL] = sort(tmp,'descend');

%Pick the best limit.
if rnkUL(2) == 1
    %nothing need be done
else
    UL=UL2;
end

% For "b2"
M2=sqrt(1+(b2^2));
% Note that the reference angle for cos(UL) = 1/M2 & sin(UL) = b2/M2

% Find the "correct" quadrant for angle denoted by acos(1/M1) & asin(b1/M1)
[LL]=circ_CorrAng(1/M2,b2/M2);

% Find the second option for the lower limits - see Zar pg 614
if LL+pi > 2*pi
    LL2=(LL+pi)-2*pi;
else
    LL2=LL+pi;
end

% Find the distances between mu2 and LL/LL2
% Backward LL:
if mu2 < LL
    bLL=((2*pi)-LL)+mu2;
elseif mu2 > LL
    bLL=abs(LL-mu2);
end
% Backward LL2:
if mu2 < LL2
    bLL2=((2*pi)-LL2)+mu2;
elseif mu2 > LL2
    bLL2=abs(LL2-mu2);
end

% Sort those distances to find the smallest 
tmp=[bLL,bLL2];
[~,rnkLL] = sort(tmp,'descend');


%Pick the best limit.
if rnkLL(2) == 1
    % nothing need be done
else
    LL=LL2;
end

% Check for imaginary confidence limits (when input data is not
% directional) and error with information for the user.
if imag(UL)~=0
    error("Calculated confidence limits are imaginary. This occurs when data is not significantly directional. Check directionality before calculating second order confidence limits.")
end
end