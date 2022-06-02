% function mustHaveSameNaN(a,b)
%
% Function validates that the two input arguments are the same length and
% have NaN values in the same positions.
%
% INPUT
% a: vector double input.
% b: vector double input.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [April 27, 2022]
%__________________________________________________________________________

function mustHaveSameNaN(a,b)

arguments
    a (1,:) double
    b (1,:) double
end

% test that NaNs are at the same position
if length(isnan(a))~=length(isnan(b))
    eidType = 'mustHaveSameNaN:notSameLength';
    msgType = 'Inputs must be the same length.';
    throwAsCaller(MException(eidType,msgType))
end

% test that NaNs are at the same position
if ~all(isnan(a)==isnan(b))
    eidType = 'mustHaveSameNaN:notSameNaN';
    msgType = 'Inputs must have NaNs in the same position.';
    throwAsCaller(MException(eidType,msgType))
end
end