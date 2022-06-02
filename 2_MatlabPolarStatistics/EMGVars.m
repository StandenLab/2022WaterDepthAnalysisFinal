%% README
% This script is written to work after data has been turned into polar
% coordinates. Thus, you're starting with all of the individual timings
% already colated into files saved for each condition.
%
% Water depth treatments are referred to by abbreviations that are viable
% matlab names as follows:
% 3.0 BD - three
% 1.1 BD - above
% 1.0 BD - at
% 0.9 BD - below
% 0.7 BD - eyes
%
% The script below runs the reported polar statistics for the timing of
% body muscle activity onset and offset.
%% Run the inital stats:
% Hermans-Rasson, Rayleigh's test, Kupier's test for a von Mises
% distribution. Also calculate angular mean, angular variance, and
% confidence limits (when appropriate).

% Body muscle activity onset and offset
clear
clc
[pn]=uigetdir("","Select your EMGPolar Directory");
flls=["_three.mat",...
    "_above.mat","_at.mat","_below.mat",...
    "_eyes.mat"];
bodycols=["e1","e2","e3","e4","e5","e6","e7"];  % electrode positions

wb = waitbar(0,'Initializing...');              % generate waitbar
for f=1:length(flls)
    waitbar(f/(length(flls)),wb,{'Checking for directional data:','Body Muscle Activity'});

    cond=extractAfter(extractBefore(flls(f),".mat"),"_");                               % pull out condition
    working=cleanload(append(pn,"\",flls(f)));                                          % load in the data
    [EPS.(cond).BonT,EPS.(cond).BonD]=WDstats(working.body.onsetbytailALL,bodycols);    % muscle activity onset
    [EPS.(cond).BoffT,EPS.(cond).BoffD]=WDstats(working.body.offsetbytailALL,bodycols); % muscle activity offset
    
    clear working cond
end
close(wb)

% Save out preliminary circular statistics. Output data is sorted by
% condition inside the structure "EPS". Inside each condition structure,
% variables are as follows:
% Bon - body muscle activity onset
% Boff - body muscle activity offset
%
% Variables with "T" at the end are output statistical tables. Variables
% with "D" at the end are the raw data. Body positions in the raw data are
% across columns.

save("WaterDepthManuscript_EMGPolarStats.mat","EPS");   % save data in current working directory.
%% Run the specific stats for each variable
clear
clc
load("WaterDepthManuscript_EMGPolarStats.mat")

%Run comparisons between means.
%Onset
[on,condon]=alldata(EPS,'BonD');                                                        %Pull together all the data in each condition
[onF,onpair,onpairB]=wwtestandmultcomp(on(1:5),condon(1:5),["e1","e2","e3","e4","e5"]); %Overall and pairwise watson williams tests
[onlettersB]=pairwiseletters(onpairB);                                                  %Letter report

%Offset
[off,condoff]=alldata(EPS,'BoffD');                                                          %Pull together all the data in each condition
[offF,offpair,offpairB]=wwtestandmultcomp(off(1:5),condoff(1:5),["e1","e2","e3","e4","e5"]); %Overall and pairwise watson williams tests
[offlettersB]=pairwiseletters(offpairB);                                                     %Letter report