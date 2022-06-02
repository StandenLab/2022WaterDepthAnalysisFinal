% function [ang]=circ_CorrAng(cosA,sinA)
%
% This function determines the quadrant postion of your input angle based
% on the cos and sin of that angle. Spits out the angle in the correct
% quadrant. Input is in RADIANS.
%
% INPUT
% cosA: cosine of your angle of interest
% sinA: sine of your angle of interest
%
% OUTPUT
% ang: your angle in the correct quadrant, based upon the sign of the
%   input cosine and sine values. Outpu in RADIANS.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [November 21, 2018]
%__________________________________________________________________________

function [ang]=circ_CorrAng(cosA,sinA)

% Determine the quadrant and output the mean angle.
if cosA > 0 && sinA > 0 %i.e. Quadrant 1
    disp('Mean angle in Q1. Angle determined by acos')
    ang=acos(cosA);
elseif cosA == 0 && sinA == 1  %i.e. pi/2
    disp('Mean angle = pi/2.')
    ang=pi/2;
elseif cosA < 0 && sinA > 0  %i.e. Quadrant 2
    disp('Mean angle in Q2. Angle determined by acos')
    ang=acos(cosA);
elseif cosA == -1 && sinA == 0  %i.e. pi
    disp('Mean angle = pi.')
    ang=pi;
elseif cosA < 0 && sinA < 0  %i.e. Quadrant 3
    disp('Mean angle in Q3. Angle determined as pi + atan')
    ang=pi+atan(sinA/cosA);
elseif cosA == 0 && sinA == -1  %i.e. 3pi/2
    disp('Mean angle = 3*pi/2.')
    ang=3*pi/2;
elseif cosA >0 && sinA < 0  %i.e. Quadrant 4
    disp('Mean angle in Q4. Angle determined as 2*pi - acos')
    ang=2*pi-acos(cosA);
end
end