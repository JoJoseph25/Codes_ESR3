% This is the main code that will call underlying functions to carry out
% various parameter calcultions and plotting.

% Clear workspace, clear command window, close any open figure/window
clearvars; clc; close all;

% DIRECTORY VARIABLE DEFINITONS
% (Main Directoy - Code Repositiory with all functions)

main_workdir = 'G:\Socrates\Codes'; % Joel Computer Working Directory

% main_workdir= 'D:\LabDrive\Lab drive new\Joel\Codes'; % Lab Computer Working Directory

cd(main_workdir); % Change to Working Directory

%Experiment Directory (where all data folder are stored)
exp_workdir='G:\Socrates\Exp_Data\Henrique_data\Exp_data'; %Joel Computer Experiment Data Location

% exp_workdir='D:\LabDrive\Lab drive new\Joel\Exp_Data\Henrique_data\Exp_data'; %Lab Computer Experiment Data Location

% SUBJECT ID & EXPERIMENTAL TRIALS

% Subject ID numbers

subject_num={...
%              '1005';'1006';'1007';'1008';...
%              '1009';'1010';'1011';'1012';...
%              '1013';'1014';'1015';'1016';...
%              '1017';'1018';'1019';'1020';...
%              '1021';'1022';'1023';
             '1024';...
             '1025';'1026';'1027';'1028';...
             };

% subject_num={'1007'};

% Experiment trail number

% exp_id ={'dy5'};

exp_id = {'dy4';'dy5';'dy6';'dy7';'dy8'};

% Exceptions

exceptions ={...
            '1020_dy5'; % COP filter check for this particular case
            '1013_dy8';'1025_dy8';'1027_dy6';
            };

% INITIALIZATION SECTION
% Initialize all frequency and file specific parameters

[init_struct_incmp] = initialize; % function call to initiate

% USER INPUTSy
[init_struct] = take_user_inp(init_struct_incmp);

%LOOPZ
for sub=1:length(subject_num)   % Main for loop going over all subjects
    for exp=1:length(exp_id)     % Second for loop going over all exp cond

        % Collect subject number from all subjects list

        x_num=subject_num{sub}; % Current Subject ID

        tr_num=exp_id{exp}; % Current Experiment Condition Number

        %% CHANGE DIRECTORY TO SUBJECT SPECIFIC FOLDER

        cd(main_workdir);

        % Function to change directory and make new folder
        [sub_file, p_num] = make_new_folder(x_num,tr_num,main_workdir,exp_workdir);

        clear x_num tr_num

        %% EXCEPTION
        % for data files that don't exist

        if ismember(p_num,exceptions)
            disp('Exception');
            continue;
        end

        %%  FETCH FILES

        % Function to fetch marker and forceplate data from files
        [FP_data, Mrk_Data,Static_Data] = fetchfiles(p_num);

        rmpath(sub_file) % remove subject file from path

        %% FIRST STEP DETECTION

        % Function to detect first step and plots
        [ind_start,ind_end] = first_step_detect(init_struct,p_num,Mrk_Data,FP_data);

        %% STANDING COP

        % Function to calculating standing COP parameters and plots
%         [struct_stand_COP] = stand_COP(init_struct,p_num,Mrk_Data,FP_data,ind_start,ind_end);
% 
%         clear struct_stand_COP

        %% STANDING ANGLES
        
        % Function to calculating standing angle parameters and plots
%         [struct_stand_angle] = stand_angle(init_struct,p_num,Mrk_Data,ind_start,ind_end);

%         clear struct_stand_angle lwrst_indx_clust rwrst_indx_clust

        %% GAIT INITIALIZATION

        % Function to calculating gait initialization parameters and plots
%         [struct_gait_init] = gait_init_calc(init_struct,p_num,ind_end,FP_data,Mrk_Data,Static_Data);
% 
%         clear struct_gait_init

        %% WALKING STEPS
        
        % Function to calculating walking back step parameters and plots
        [struct_walk_back_steps] = walk_back_steps(init_struct,p_num,ind_end,Mrk_Data,Static_Data);
        
        clear struct_walk_back_steps
        
        %% WALKING ANGLES
        
        % Function to calculating walking back step parameters and plots
        [struct_walk_back_angle] = walk_back_angle(init_struct,p_num,ind_end,Mrk_Data);
        
        clear struct_walk_back_angle

        %% STRUCT COMBINE
        
        % Function to save all the calculated parameters in one excel/mat file
        comb_struct = struct_combine(init_struct,p_num);
       
        %% Clear all unnecessary variables

        clearvars -except main_workdir exp_workdir subject_num exp_id ...
                          init_struct sub sub_file exp exceptions;

    end % end of secondary for loop
    
    % Function to save all trial data in one excel file
    Subj_Table = subj_file_save(sub_file,subject_num{sub});
    
    clear Subj_Table
    
end % end of main loop

%Function to save all exp data in one excel file
Exp_Table = exp_file_save(exp_workdir,subject_num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%