function [sub_file,p_num] = make_new_folder(x_num,tr_num,main_workdir,exp_workdir)
% function [FP_data, Mrk_Data] = fetchfiles(p_num)
% 
%  This function changes folder and makes new folder based on 
%  subject ID + Experiment trial.
% 
%   INPUT:  x_num - subject ID (string) 
%           tr_num - Experiment trial number (string)
%           main_workdir - Main working directory (string)
%           exp_workdir - Expeiment directory based on subject ID and exp
%                         trial no (string)
%
%   OUTPUT: sub_file - Subject file path (string)
%           p_num - Subject ID + Experiment condition (string)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% CREATE SUJECT FILE

% Create subject file name to navigate to directory
sub_file=strcat(exp_workdir,'\',x_num);

% Subject ID + Experiment Condition
p_num=strcat(x_num,'_',tr_num);

cd(sub_file) % change directory to subject specific file  

%% ADD TO PATH

addpath(main_workdir)   % add main directory to path
addpath(sub_file)   % add subject-specific file to path

%% CREATE NEW FOLDER
if ~exist(p_num, 'dir') % Check if folder exist
    mkdir(p_num); % make new folder
    disp("Made Folder"); % display if folder made
end

cd(p_num); % Change directory to new folder

%% DISPLAY

disp(['Filename: ',p_num,' created']); % Display filename to command line

end