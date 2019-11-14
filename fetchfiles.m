function [FP_data, Mrk_Data] = fetchfiles(p_num)
% function [FP_data, Mrk_Data] = fetchfiles(p_num)
% 
%  This function imports marker and force plate data based on
%  the subject id and experiment trial number.
% 
%   INPUT:  p_num - String cotaining subject ID 
%                   and Experiment trial number (String)
%     
%   OUTPUT: FP_data - Force plate data (table)
%           Mrk_Data - Marker Data (table)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%%  FILENAME CREATION

filename_mrk = strcat(p_num,'.tsv'); % Force Plate Data Filename
filename_FP = strcat(p_num,'_f_1.tsv'); % Marker Data Filename 

%% IMPORT FORCE PLATE DATA

opts_FP = detectImportOptions(filename_FP,'FileType','text'); % Force plate data import option definition
FP_data = readtable(filename_FP,opts_FP); % Force plate data

%% IMPORT FORCE PLATE DATA

opts_mrk = detectImportOptions(filename_mrk,'FileType','text'); % marker data import option definition
opts_mrk.VariableNamesLine=11; % Where the variable names are located
opts_mrk.DataLine=12; % Line number where the table is located
Mrk_Data = readtable(filename_mrk,opts_mrk,'ReadVariableNames',true); % Marker Data 

%% DISPLAY

disp('Files Fetched'); % print out to command line

end