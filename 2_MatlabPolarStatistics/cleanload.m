% function [output]=cleanload(filename)
%
% This function loads in a structure from a filename to a new variable name
% WITHOUT nesting the saved structure name in the file inside that new name.
% This is useful for iterating across files when each file has a similarly
% nested structure, but with a different top level name.
%
% INPUT
% filename: the name of the file to load.
%
% OUTPUT
% output: the structure output.
%
%__________________________________________________________________________
% Written by: Keegan Lutek [December 16, 2021]
%__________________________________________________________________________

function [output]=cleanload(filename)
datastruct = load(filename); %Load file
storedvars = fieldnames(datastruct); %Get the fieldnames
thisdepth = storedvars{1}; %find the first fieldname
output = datastruct.(thisdepth); %place into output
end