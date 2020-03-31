function [FP_data, Mrk_Data,Static_Data] = fetchfiles(p_num)
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
%           Static_Data - Static File parameters (struct)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%%  FILENAME CREATION

filename_mrk = strcat(p_num,'.tsv'); % Force Plate Data Filename
filename_FP = strcat(p_num,'_f_1.tsv'); % Marker Data Filename 

p_num_array = strsplit(p_num,'_');
sub_id = p_num_array{1};
filename_static = strcat(sub_id,'_static.tsv'); % Static Filename

%% IMPORT FORCE PLATE DATA

opts_FP = detectImportOptions(filename_FP,'FileType','text'); % Force plate data import option definition
FP_data = readtable(filename_FP,opts_FP); % Force plate data

%% IMPORT MARKER DATA

opts_mrk = detectImportOptions(filename_mrk,'FileType','text'); % marker data import option definition
opts_mrk.VariableNamesLine=11; % Where the variable names are located
opts_mrk.DataLine=12; % Line number where the table is located
Mrk_Data = readtable(filename_mrk,opts_mrk,'ReadVariableNames',true); % Marker Data 

%% IMPORT STATIC FILE DATA
opts_static = detectImportOptions(filename_static,'FileType','text'); % marker data import option definition
opts_static.VariableNamesLine=11; % Where the variable names are located
opts_static.DataLine=12; % Line number where the table is located
Static_Tbl = readtable(filename_static,opts_static,'ReadVariableNames',true); % Marker Data 

% Structure
Static_Data =struct;

rleg_len = (Static_Tbl.RHIPZ-Static_Tbl.RANKZ);
lleg_len = (Static_Tbl.LHIPZ-Static_Tbl.LANKZ);
Static_Data.leg_len = mean(mean([rleg_len,lleg_len],2))/1000;

Static_Data.step_width = mean(sqrt((Static_Tbl.LANKY-Static_Tbl.RANKY).^2))/1000;

XCordinates =[Static_Tbl.RHEELX,Static_Tbl.RANKX,Static_Tbl.RTOEX,Static_Tbl.LTOEX,Static_Tbl.LANKX,Static_Tbl.LHEELX,Static_Tbl.RHEELX];
YCordinates =[Static_Tbl.RHEELY,Static_Tbl.RANKX,Static_Tbl.RTOEY,Static_Tbl.LTOEY,Static_Tbl.LANKY,Static_Tbl.LHEELY,Static_Tbl.RHEELY];
Static_Data.BoS = mean(polyarea(XCordinates',YCordinates'),2)/1000000;

Static_Data.Leg_norm = Static_Data.leg_len;
Static_Data.Time_norm = sqrt(Static_Data.leg_len/9.81);
Static_Data.Vel_norm = sqrt(Static_Data.leg_len*9.81);

%% DISPLAY

disp('Files Fetched'); % print out to command line

%% CHECK IF MISSING MARKERS
nan_vec = find(Mrk_Data{:,:} == 0);

if isempty(nan_vec)
    fill_data_mat = Mrk_Data{:,3:end};
else
    fill_data_mat = fill_data(p_num);
end

%% REPLACE GAPS
fill_data_mat(isnan(fill_data_mat)) = 0;

Mrk_Data{:,3:end} = fill_data_mat;

end