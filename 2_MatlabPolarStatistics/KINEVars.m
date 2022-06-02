%% README
% This script is written to work after data has been turned into polar
% coordinates. Thus, you're starting with all of the individual timings
% already colated into files saved for each condition.
%
% Water depth treatments are referred to by abbreviations that are viable
% matlab names as follows:
% 3.0 BD - three
% 2.0 BD - two
% 1.1 BD - above
% 1.0 BD - at
% 0.9 BD - below
% 0.7 BD - eyes
% 0.1 BD - mouth
% 0.0 BD - dry
%
% The script below runs the reported polar statistics for the timing of
% body undulation, pectoral fin cycle and nose elevation.
%% Run the inital stats: 
% Hermans-Rasson, Rayleigh's test, Kupier's test for a von Mises
% distribution. Also calculate angular mean, angular variance, and
% confidence limits (when appropriate).

% Body Amplitude Timing and Fin Start Cycle
clear
clc
[pn]=uigetdir("","Select KINEPolar folder");
% These are the files with the polar timings for each variable
flls=["three.mat","two.mat",...
    "above.mat","at.mat","below.mat",...
    "eyes.mat","mouth.mat",...
    "dry.mat"];
bodycols=["H","20","40","60","80","T"]; % Body positions
fincols=["L","R"];                      % Fins

wb = waitbar(0,'Initializing...');      % Generate waitbar
for f=1:length(flls)
    % Update waitbar
    waitbar(f/(length(flls)),wb,{'Checking for directional data:','Body undulation and fin movement'});
    
    % Load the appropriate file
    cond=extractBefore(flls(f),".mat");                                                     % Pull out condition
    working=cleanload(append(pn,"\",flls(f)));                                              % Load in the data
    
    % Run preliminary statistics to check for von Mises distribution, run
    % Rayleigh's or Hermans-Rasson tests as appropriate and caluculate the
    % angular mean and variance as appropriate.
    [KPS.(cond).BstartT,KPS.(cond).BstartD]=WDstats(working.body.startbytailALL,bodycols);  % Body undulation maximum left by tail
    [KPS.(cond).BmidT,KPS.(cond).BmidD]=WDstats(working.body.midbytailALL,bodycols);        % Body undulation maximum right by tail
    [KPS.(cond).FstartT,KPS.(cond).FstartD]=WDstats(working.fin.startbyfinALL,fincols);     % Pectoral fin start of cycle by right pectoral fin
    clear working cond
end
close(wb)

% Nose Elevation
clearvars -except KPS
clc
[pn]=uigetdir("","Select ELAVPolar folder");
% These are the files with the polar timings for nose elevation
flls=["above.mat","at.mat","below.mat",...
    "eyes.mat","mouth.mat",...
    "dry.mat"];

wb = waitbar(0,'Initializing...');      % Generate waitbar
for f=1:length(flls)
    % Update waitbar
    waitbar(f/(length(flls)),wb,{'Checking for directional data:','Nose elevation'});
    
    % Load the appropriate file
    cond=extractBefore(flls(f),".mat");                                         % Pull out the name of the condition based on the filename.
    working=cleanload(append(pn,"\",flls(f)));                                  % Load in the data
    
    % Run preliminary statistics to check for von Mises distribution, run
    % Rayleigh's or Hermans-Rasson tests as appropriate and caluculate the
    % angular mean and variance as appropriate.
    [KPS.(cond).LLT,KPS.(cond).LLD]=WDstats(working.body.leftbyleftALL,"L");    % Nose elevation to the left by left pectoral fin stroke
    [KPS.(cond).RRT,KPS.(cond).RRD]=WDstats(working.body.rightbyrightALL,"R");  % Nose elevation to the right by right pectoral fin stroke
    clear working cond
end
close(wb)

% Save out preliminary circular statistics. Output data is sorted by
% condition inside the structure "KPS". Inside each condition structure,
% variables are as follows:
% Bstart - body maximum amplitude to the left
% Bmid - body maximum amplitude to the right
% Fstart - start of the pectoral fin cycle
% LL - nose maximum elevation to the left
% RR - nose maximum elevation to the right
%
% Variables with "T" at the end are output statistical tables. Variables
% with "D" at the end are the raw data. Body and fin positions in the raw
% data are across columns as specified in 'bodycols' and 'fincols'.

save("WaterDepthManuscript_KinePolarStats.mat","KPS");  % save data in current working directory
%% Run the specific tests for each variable
clear
clc
load("WaterDepthManuscript_KinePolarStats.mat")

% Pectoral fins don't need anything further run.

% Body needs to run a comparison between means.
% Maximum amplitude to the left
[start,condstart]=alldata(KPS,'BstartD');                                                       % Pull together all the data in each condition
[startF,startpair,startpairB]=wwtestandmultcomp(start,condstart,["H","20","40","60","80","T"]); % Overall and pairwise watson williams tests
[startlettersB]=pairwiseletters(startpairB);                                                    % Letter report
% Maximum amplitude to the right
[mid,condmid]=alldata(KPS,'BmidD');                                                             % Pull together all the data in each condition
[midF,midpair,midpairB]=wwtestandmultcomp(mid,condmid,["H","20","40","60","80","T"]);           % Overall and pairwise watson williams tests
[midlettersB]=pairwiseletters(midpairB);                                                        % Letter report

% Nose elevation needs to run a comparison between angular dispersion.
% Nose elevation to the left
[left,condleft,leftdif]=alldata(KPS,'LLD','angdif',1);                                      % Pull together all the data in each condition
[FpL,FtblL,Lstats]=kruskalwallis(leftdif{1},condleft{1});                                   % Test for differences in angular variance
LmultcompE = splitvars(table(multcompare(Lstats)));                                         % Post-hoc multiple comparisons
LmultcompE.Properties.VariableNames={'group1','group2','UL','diff','LL','p-value'};         % Rename the table columns
LmultcompE.("p-valueB")=LmultcompE.("p-value")*size(LmultcompE,1);                          % Calculate the bonferroni corrected p-values
LmultcompE.("p-valueB")(LmultcompE.("p-valueB")>1)=1;                                       % Replace anything larger than 1 with 1 cause probabilities above 1 make no sense
[leftlettersE]=letterreport(LmultcompE.("p-valueB"),LmultcompE.group1,LmultcompE.group2);   % Letter report
% Nose elevation to the right
[right,condright,rightdif]=alldata(KPS,'RRD','angdif',1);                                   % Pull together all the data in each condition
[FpR,FtblR,Rstats]=kruskalwallis(rightdif{1},condright{1});                                 % Test for differences in angular variance
RmultcompE = splitvars(table(multcompare(Rstats)));                                         % Post-hoc multiple comparisons
RmultcompE.Properties.VariableNames={'group1','group2','UL','diff','LL','p-value'};         % Rename the table columns
RmultcompE.("p-valueB")=RmultcompE.("p-value")*size(RmultcompE,1);                          % Calculate the bonferroni corrected p-values
RmultcompE.("p-valueB")(RmultcompE.("p-valueB")>1)=1;                                       % Replace anything larger than 1 with 1 cause probabilities above 1 make no sense
[rightlettersE]=letterreport(RmultcompE.("p-valueB"),RmultcompE.group1,RmultcompE.group2);  % Letter report