function [Fill_Data] = fill_data(p_num)
% function function [Fill_Data] = fill_data(filename)
% 
%  This function imports marker data based on the 
%  filename and fills it in based on Gloersen et al. 2016.
% 
%   INPUT:  p_num - String cotaining subject ID 
%                   and Experiment trial number (String)
%     
%   OUTPUT: Fill_Data - Filled Marker data (matrix)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% EXCEPTIONS

exceptions ={'1006_dy5';'1011_dy7';'1012_dy7';'1016_dy5';'1028_dy8';};

if ismember(p_num,exceptions)
    filename = strcat(p_num,'.c3d');
else
    filename = strcat(p_num,'.tsv');
end

%% Initialize and load data

disp("Filling Missing Markers: ");

dataFile = filename;
connectionFile = 'connect_file.txt';

%Load data (requires the MoCap toolbox)
original = mcread(dataFile);

%Remove invalid markers in hdm files
markernames = original.markerName;
ids = find(~ismember(markernames,{'*0','*1','*2'}));
original = mcgetmarker(original,ids);

%% Original data

%MoCap display parameters
p = mcinitanimpar;
p = mccreateconnmatrix(connectionFile,p);

incomplete = original;

%% Recovery examples

% PMA of different individual models

%Different combinations may be tried to extract the best results

%options.saveastsv = 0;%save recovered sequence as tsv file
%options.recursivefilling = 0;%use or not recusrivefilling
options.method1 = 1;%Local interpolation
options.method2 = 1;% Local polynomial regression
options.method3 = 1;%Local GRNN
options.method4 = 1;%Global weighted linear regression
options.method5 = 1;%Gloersen et al. 2016
options.advancedordering = 1;
options.spaceconstraint = 1;%use or not spaceconstraint
options.timeconstraint = 1;%use or not timeconstraint
options.filtering = 1;%use or not timeconstraint
options.quiet = 0;%avoid console output
options.presenceMin = 30;%threshold (in % of available frames) under which discard some markers

recovered = mcrecovery(incomplete,options);

%% OUTPUT DATA FILE

Fill_Data = recovered.data;

end