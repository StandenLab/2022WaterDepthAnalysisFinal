% function [pval,teststat,critval,norm] = circ_rayleigh2(data,r)
%
% This function determines whether a set of second order means is
% direction. First checks if rectangular coordinates are bivarate normally
% distributed; if so, calculates F statistic to check for directionality of
% set of mean angles. If not, calculates R statistic based on non-parametric
% test for directionality. From Zar (1999) pgs 638-640.
%
% Ho: There is no mean population direction (uniform distribution of the
% values around the circle).
%
% INPUT
% data: set of mean angles. Data in a single row.
% r: set of average vector lengths for your set of mean angles.
%
% OUTPUT
% pval: p-value for the test of directionality against Ho.
% teststat: test statistic for test of directionality.
% critval: critical value against which the test statistic is compared.
% norm: are the rectangular coordinates for your mean BOTH normally
%   distributed. 'Normal Data' means you identified both X and Y as normal.
%   'Non-Normal Data' means you identified at least one of X and Y as
%   non-normal.
% 
% REQUIRES
% mustHaveSameNaN
% moorelookup
%
% REFERENCES
% Zar, J.H. (1999). Biostatistical Analysis. Fourth edition. Upper Saddle
%   River, New Jersey: Prentice-Hall Inc.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [November 22, 2018]
%
% Edited:
% [20220428 | KL] - cleaning function and comments
%__________________________________________________________________________

function [pval,teststat,critval,norm] = circ_rayleigh2(data,r)

arguments
    data (1,:) double {mustBeVector}
    r (1,:) double {mustBeVector, mustHaveSameNaN(data,r)}
end

%% Initial data [k,X,Y,XY]
% Set # of samples (k)
k=size(data,2);

% Calculate rectangular coordinates for each mean angle.
X=nan(1,length(data));      % preset for speed
Y=nan(1,length(data));      % preset for speed
for i=1:size(data,2)
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
%% Check for normality of X and Y.
% Compute Shapiro-Wilks Test p-value and test statistic for X and Y
[HX,pX,statX]=swtest(X);
[HY,pY,statY]=swtest(Y);

% Create probablility density curves for theoretical X and Y (based on
% their mean and standard deviation)
% Histogram plot for X
figure(1);
subplot(2,2,1)
hX=histfit(X);hold on;
hX(1).EdgeColor=[0,0,0];
hX(1).FaceColor=[0,0,0];
if HX == 1
    title(sprintf('Histogram - X \n p = %.4f; W = %.4f \n Reject H0. Data not normal',pX,statX))
else
    title(sprintf('Histogram - X \n p = %.4f; W = %.4f \n Accept H0. Data are normal',pX,statX))
end
% Histogram plot for Y
subplot(2,2,2)
hY=histfit(Y);hold on;
hY(1).EdgeColor=[0,0,0];
hY(1).FaceColor=[0,0,0];
if HY == 1
    title(sprintf('Histogram - Y \n p = %.4f; W = %.4f \n Reject H0. Data not normal',pY,statY))
else
    title(sprintf('Histogram - Y \n p = %.4f; W = %.4f \n Accept H0. Data are normal',pY,statY))
end
% QQ norm plot for X
qqX=subplot(2,2,3);
qqplot(X);
qqX.XLabel.String='';
qqX.YLabel.String='';
qqX.Title.String='QQNorm - X';
% QQ norm plot for Y
qqY=subplot(2,2,4);
qqplot(Y);
qqY.XLabel.String='';
qqY.YLabel.String='';
qqY.Title.String='QQNorm - Y';

% User makes final decision of how to treat the data.
answer = questdlg('Are BOTH X and Y normal?',...
    'Normality Check',...
    'Yes','No','Yes');

close(figure(1))
%% Test for directionality
switch answer
    case 'Yes' %i.e. BOTH X and Y are normal.
        F=((k*(k-2))/2)*(((Xbar^2*ssY)-(2*Xbar*Ybar*ssXY)+(Ybar^2*ssX))/((ssX*ssY)-ssXY^2));    % F test statistic
        Fcrit=finv(0.95,2,k-2);                                                                 % F critical value
        Fp=1-fcdf(F,2,k-2);                                                                     % F p-value
        pval=Fp;
        teststat=F;
        critval=Fcrit;
        norm='Normal Data';
    case 'No' %i.e. at least one of X and Y is non-normal.
        % Rank data.
        [~,rank] = sort(r,'ascend');
        % Compute X/Y numerators 
        numX=nan(1,length(data));      % preset for speed
        numY=nan(1,length(data));      % preset for speed
        for i=1:size(X,2)
        numX(i)=rank(i)*cos(data(i));
        numY(i)=rank(i)*sin(data(i));
        end
        % Compute nonparametric X and Y
        Xnp=sum(numX)/k;
        Ynp=sum(numY)/k;
        % Compute nonparametric test stat R'
        R=sqrt((Xnp^2+Ynp^2)/k);
        % Check where R' falls.
        [p2side, Rcrit] = moorelookup(k, R);
        pval=p2side;
        teststat=R;
        critval=Rcrit;
        norm='Non-Normal Data';
end
end