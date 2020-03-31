function [comb_struct] = struct_combine(init_struct,p_num)
% function [comb_struct] = struct_combine(p_num)
% 
%  This function combines all struct to one single structure for that
%  session.
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Number (string)
%     
%   OUTPUT: comb_struct - Session parameters (struct)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il) 

%% PATH & USER INFO

%Emotion Directory (where all emotion specific files are stored)
emo_workdir='../../../OpenSeasme'; %Joel Computer Emotion Data Location

addpath(emo_workdir);

p_num_array = strsplit(p_num,'_');

sub_id = p_num_array{1};
filename_emo = strcat('subject-',sub_id,'.csv');

trial_num = p_num_array{2};
exp_num = str2num(trial_num(3));

%% FETCH EMOTION DATA

opts_emo = detectImportOptions(filename_emo,'FileType','text'); % Emotion file import options 
Emo_data = readtable(filename_emo,opts_emo); % Emotion Table

diff = length(Emo_data.Emotion)-8; % difference between rows and exp_num
exp_num = exp_num + diff;

emo_name = string(cell2mat(Emo_data.Emotion(exp_num)));
emo_num = Emo_data.emo_code(exp_num);

rmpath(emo_workdir);

disp(strcat('Emotion-',emo_name,' Structure Combined: '));

%% LOAD & SAVE .MAT

comb_struct = struct;

comb_struct.subject_id = string(sub_id);
comb_struct.trial_num = string(trial_num);
comb_struct.emo_num =str2double(emo_num(:));
comb_struct.emo_name = emo_name;

% Stand COP
addpath('./Stand COP');

load(strcat(p_num,'_stand_COP.mat'));

comb_struct.COP_transitions = struct_stand_COP.transitions;

comb_struct.COP_full_mean_COP_X = struct_stand_COP.full_mean_COP_X;
comb_struct.COP_full_mean_COP_Y = struct_stand_COP.full_mean_COP_Y;
comb_struct.COP_full_mean_COP_XY = struct_stand_COP.full_mean_COP_XY;
comb_struct.COP_full_mean_COP_XY_ratio = struct_stand_COP.full_mean_COP_XY_ratio;

comb_struct.COP_full_std_COP_X = struct_stand_COP.full_std_COP_X;
comb_struct.COP_full_std_COP_Y = struct_stand_COP.full_std_COP_Y;
comb_struct.COP_full_std_COP_XY = struct_stand_COP.full_std_COP_XY;
comb_struct.COP_full_std_COP_XY_ratio = struct_stand_COP.full_std_COP_XY_ratio;

comb_struct.COP_full_vel_mean_COP_X = struct_stand_COP.full_vel_mean_COP_X;
comb_struct.COP_full_vel_mean_COP_Y = struct_stand_COP.full_vel_mean_COP_Y;
comb_struct.COP_full_vel_mean_COP_XY = struct_stand_COP.full_vel_mean_COP_XY;
comb_struct.COP_full_vel_mean_COP_XY_ratio = struct_stand_COP.full_vel_mean_COP_XY_ratio;

comb_struct.COP_full_vel_std_COP_X = struct_stand_COP.full_vel_std_COP_X;
comb_struct.COP_full_vel_std_COP_Y = struct_stand_COP.full_vel_std_COP_Y;
comb_struct.COP_full_vel_std_COP_XY = struct_stand_COP.full_vel_std_COP_XY;
comb_struct.COP_full_vel_std_COP_XY_ratio = struct_stand_COP.full_vel_std_COP_XY_ratio;

comb_struct.COP_full_mean_step_width = struct_stand_COP.full_mean_step_width;
comb_struct.COP_full_mean_BoS = struct_stand_COP.full_mean_BoS;

comb_struct.COP_full_std_step_width = struct_stand_COP.full_std_step_width;
comb_struct.COP_full_std_BoS = struct_stand_COP.full_std_BoS;

comb_struct.COP_clust_rel_mean_COP_X = struct_stand_COP.clust_rel_mean_COP_X;
comb_struct.COP_clust_rel_mean_COP_Y = struct_stand_COP.clust_rel_mean_COP_Y;
comb_struct.COP_clust_rel_mean_COP_XY = struct_stand_COP.clust_rel_mean_COP_XY;
comb_struct.COP_clust_rel_mean_COP_XY_ratio = struct_stand_COP.clust_rel_mean_COP_XY_ratio;

comb_struct.COP_clust_rel_std_COP_X = struct_stand_COP.clust_rel_std_COP_X;
comb_struct.COP_clust_rel_std_COP_Y = struct_stand_COP.clust_rel_std_COP_Y;
comb_struct.COP_clust_rel_std_COP_XY = struct_stand_COP.clust_rel_std_COP_XY;
comb_struct.COP_clust_rel_std_COP_XY_ratio = struct_stand_COP.clust_rel_std_COP_XY_ratio;

comb_struct.COP_clust_rel_vel_mean_COP_X = struct_stand_COP.clust_rel_vel_mean_COP_X;
comb_struct.COP_clust_rel_vel_mean_COP_Y = struct_stand_COP.clust_rel_vel_mean_COP_Y;
comb_struct.COP_clust_rel_vel_mean_COP_XY = struct_stand_COP.clust_rel_vel_mean_COP_XY;
comb_struct.COP_clust_rel_vel_mean_COP_XY_ratio = struct_stand_COP.clust_rel_vel_mean_COP_XY_ratio;

comb_struct.COP_clust_rel_vel_std_COP_X = struct_stand_COP.clust_rel_vel_std_COP_X;
comb_struct.COP_clust_rel_vel_std_COP_Y = struct_stand_COP.clust_rel_vel_std_COP_Y;
comb_struct.COP_clust_rel_vel_std_COP_XY = struct_stand_COP.clust_rel_vel_std_COP_XY;
comb_struct.COP_clust_rel_vel_std_COP_XY_ratio = struct_stand_COP.clust_rel_vel_std_COP_XY_ratio;

comb_struct.COP_clust_rel_mean_step_width = struct_stand_COP.clust_rel_mean_step_width;
comb_struct.COP_clust_rel_mean_BoS = struct_stand_COP.clust_rel_mean_BoS;

comb_struct.COP_clust_rel_std_step_width = struct_stand_COP.clust_rel_std_step_width;
comb_struct.COP_clust_rel_std_BoS = struct_stand_COP.clust_rel_std_BoS;

clear struct_stand_COP
rmpath('./Stand COP');

disp('    Stand COP Saved');

% Stand Angle
addpath('./Stand Angle');

load(strcat(p_num,'_stand_angle.mat'));

comb_struct.standing_full_mean_SHO_XY_ang = struct_stand_angle.full_mean_SHO_xy_ang;
comb_struct.standing_full_mean_SHO_XZ_ang = struct_stand_angle.full_mean_SHO_xz_ang;
comb_struct.standing_full_mean_head_ang = struct_stand_angle.full_mean_head_ang;
comb_struct.standing_full_mean_back_ang = struct_stand_angle.full_mean_back_ang;
comb_struct.standing_full_mean_rwrst_dist = struct_stand_angle.full_mean_rwrst_dist;
comb_struct.standing_full_mean_lwrst_dist = struct_stand_angle.full_mean_lwrst_dist;

comb_struct.standing_full_std_SHO_XY_ang = struct_stand_angle.full_std_SHO_xy_ang;
comb_struct.standing_full_std_SHO_XZ_ang = struct_stand_angle.full_std_SHO_xz_ang;
comb_struct.standing_full_std_head_ang = struct_stand_angle.full_std_head_ang;
comb_struct.standing_full_std_back_ang = struct_stand_angle.full_std_back_ang;
comb_struct.standing_full_std_rwrst_dist = struct_stand_angle.full_std_rwrst_dist;
comb_struct.standing_full_std_lwrst_dist = struct_stand_angle.full_std_lwrst_dist;

comb_struct.standing_clust_rel_mean_SHO_XY_ang = struct_stand_angle.clust_rel_mean_SHO_xy_ang;
comb_struct.standing_clust_rel_mean_SHO_XZ_ang = struct_stand_angle.clust_rel_mean_SHO_xz_ang;
comb_struct.standing_clust_rel_mean_head_ang = struct_stand_angle.clust_rel_mean_head_ang;
comb_struct.standing_clust_rel_mean_back_ang = struct_stand_angle.clust_rel_mean_back_ang;
comb_struct.standing_clust_rel_mean_rwrst_dist = struct_stand_angle.clust_rel_mean_rwrst_dist;
comb_struct.standing_clust_rel_mean_lwrst_dist = struct_stand_angle.clust_rel_mean_lwrst_dist;

comb_struct.standing_clust_rel_std_SHO_XY_ang = struct_stand_angle.clust_rel_std_SHO_xy_ang;
comb_struct.standing_clust_rel_std_SHO_XZ_ang = struct_stand_angle.clust_rel_std_SHO_xz_ang;
comb_struct.standing_clust_rel_std_head_ang = struct_stand_angle.clust_rel_std_head_ang;
comb_struct.standing_clust_rel_std_back_ang = struct_stand_angle.clust_rel_std_back_ang;
comb_struct.standing_clust_rel_std_rwrst_dist = struct_stand_angle.clust_rel_std_rwrst_dist;
comb_struct.standing_clust_rel_std_lwrst_dist = struct_stand_angle.clust_rel_std_lwrst_dist;

clear struct_stand_angle
rmpath('./Stand Angle');

disp('    Stand Angles Saved');

% Gait Initialization
addpath('./Gait Init');

load(strcat(p_num,'_gait_init.mat'));

comb_struct.gait_init_foot = struct_gait_init.foot_gait_init;

comb_struct.gait_init_swing_len = struct_gait_init.swing_len_gait_init;
comb_struct.gait_init_swing_duration = struct_gait_init.swing_duration_gait_init;
comb_struct.gait_init_swing_speed = struct_gait_init.swing_speed_gait_init;
comb_struct.gait_init_swing_height = struct_gait_init.swing_height_gait_init;
comb_struct.gait_init_swing_width = struct_gait_init.swing_width_gait_init;

comb_struct.first_step_foot = struct_gait_init.first_step;

comb_struct.first_step_swing_len = struct_gait_init.swing_len_first_step;
comb_struct.first_step_swing_duration = struct_gait_init.swing_duration_first_step;
comb_struct.first_step_swing_speed = struct_gait_init.swing_speed_first_step;
comb_struct.first_step_swing_height = struct_gait_init.swing_height_first_step;
comb_struct.first_step_swing_width = struct_gait_init.swing_width_first_step;

comb_struct.gait_init_norm_swing_len = struct_gait_init.norm_swing_len_gait_init;
comb_struct.gait_init_norm_swing_duration = struct_gait_init.norm_swing_duration_gait_init;
comb_struct.gait_init_norm_swing_speed = struct_gait_init.norm_swing_speed_gait_init;

comb_struct.first_step_norm_swing_len = struct_gait_init.norm_swing_len_first_step;
comb_struct.first_step_norm_swing_duration = struct_gait_init.norm_swing_duration_first_step;
comb_struct.first_step_norm_swing_speed = struct_gait_init.norm_swing_speed_first_step;

clear struct_gait_init
rmpath('./Gait Init');

disp('    Gait Initilization Saved');

% Walking Steps
addpath('./Walking Steps');

load(strcat(p_num,'_walk_back_steps.mat'));

comb_struct.walking_right_step_count = struct_walk_back_steps.right_step_count;

comb_struct.walking_mean_right_step_len = struct_walk_back_steps.right_step_len_mean;
comb_struct.walking_right_mean_step_duration = struct_walk_back_steps.right_step_duration_mean;
comb_struct.walking_right_mean_step_speed = struct_walk_back_steps.right_step_speed_mean;
comb_struct.walking_right_mean_step_height = struct_walk_back_steps.right_step_height_mean;
comb_struct.walking_right_mean_step_width = struct_walk_back_steps.right_step_width_mean;

comb_struct.walking_right_mean_swing_len = struct_walk_back_steps.right_swing_len_mean;
comb_struct.walking_right_mean_swing_duration = struct_walk_back_steps.right_swing_duration_mean;
comb_struct.walking_right_mean_swing_speed = struct_walk_back_steps.right_swing_speed_mean;
comb_struct.walking_right_mean_swing_ratio = struct_walk_back_steps.right_swing_ratio_mean;

comb_struct.walking_right_mean_stance_len = struct_walk_back_steps.right_stance_len_mean;
comb_struct.walking_right_mean_stance_duration = struct_walk_back_steps.right_stance_duration_mean;
comb_struct.walking_right_mean_stance_speed = struct_walk_back_steps.right_stance_speed_mean;
comb_struct.walking_right_mean_stance_ratio = struct_walk_back_steps.right_stance_ratio_mean;

comb_struct.walking_right_mean_temp_ratio = struct_walk_back_steps.right_temp_ratio_mean;

comb_struct.walking_right_step_count = struct_walk_back_steps.right_step_count;

comb_struct.walking_right_std_step_len = struct_walk_back_steps.right_step_len_std;
comb_struct.walking_right_std_step_duration = struct_walk_back_steps.right_step_duration_std;
comb_struct.walking_right_std_step_speed = struct_walk_back_steps.right_step_speed_std;
comb_struct.walking_right_std_step_height = struct_walk_back_steps.right_step_height_std;
comb_struct.walking_right_std_step_width = struct_walk_back_steps.right_step_width_std;

comb_struct.walking_right_std_swing_len = struct_walk_back_steps.right_swing_len_std;
comb_struct.walking_right_std_swing_duration = struct_walk_back_steps.right_swing_duration_std;
comb_struct.walking_right_std_swing_speed = struct_walk_back_steps.right_swing_speed_std;
comb_struct.walking_right_std_swing_ratio = struct_walk_back_steps.right_swing_ratio_std;

comb_struct.walking_std_right_stance_len = struct_walk_back_steps.right_stance_len_std;
comb_struct.walking_right_std_stance_duration = struct_walk_back_steps.right_stance_duration_std;
comb_struct.walking_right_std_stance_speed = struct_walk_back_steps.right_stance_speed_std;
comb_struct.walking_right_std_stance_ratio = struct_walk_back_steps.right_stance_ratio_std;

comb_struct.walking_right_std_temp_ratio = struct_walk_back_steps.right_temp_ratio_std;


comb_struct.walking_norm_right_mean_step_len = struct_walk_back_steps.norm_right_step_len_mean;
comb_struct.walking_norm_right_mean_step_duration = struct_walk_back_steps.norm_right_step_duration_mean;
comb_struct.walking_norm_right_mean_step_speed = struct_walk_back_steps.norm_right_step_speed_mean;
comb_struct.walking_norm_right_mean_swing_len = struct_walk_back_steps.norm_right_swing_len_mean;
comb_struct.walking_norm_right_mean_swing_duration = struct_walk_back_steps.norm_right_swing_duration_mean;
comb_struct.walking_norm_right_mean_swing_speed = struct_walk_back_steps.norm_right_swing_speed_mean;
comb_struct.walking_norm_right_mean_stance_len = struct_walk_back_steps.norm_right_stance_len_mean;
comb_struct.walking_norm_right_mean_stance_duration = struct_walk_back_steps.norm_right_stance_duration_mean;
comb_struct.walking_norm_right_mean_stance_speed = struct_walk_back_steps.norm_right_stance_speed_mean;

comb_struct.walking_norm_right_std_step_len = struct_walk_back_steps.norm_right_step_len_std;
comb_struct.walking_norm_right_std_step_duration = struct_walk_back_steps.norm_right_step_duration_std;
comb_struct.walking_norm_right_std_step_speed = struct_walk_back_steps.norm_right_step_speed_std;
comb_struct.walking_norm_right_std_swing_len = struct_walk_back_steps.norm_right_swing_len_std;
comb_struct.walking_norm_right_std_swing_duration = struct_walk_back_steps.norm_right_swing_duration_std;
comb_struct.walking_norm_right_std_swing_speed = struct_walk_back_steps.norm_right_swing_speed_std;
comb_struct.walking_norm_right_std_stance_len = struct_walk_back_steps.norm_right_stance_len_std;
comb_struct.walking_norm_right_std_stance_duration = struct_walk_back_steps.norm_right_stance_duration_std;
comb_struct.walking_norm_right_std_stance_speed = struct_walk_back_steps.norm_right_stance_speed_std;



comb_struct.walking_left_step_count = struct_walk_back_steps.left_step_count;

comb_struct.walking_left_mean_step_len = struct_walk_back_steps.left_step_len_mean;
comb_struct.walking_left_mean_step_duration = struct_walk_back_steps.left_step_duration_mean;
comb_struct.walking_left_mean_step_speed = struct_walk_back_steps.left_step_speed_mean;
comb_struct.walking_left_mean_step_height = struct_walk_back_steps.left_step_height_mean;
comb_struct.walking_left_mean_step_width = struct_walk_back_steps.left_step_width_mean;

comb_struct.walking_left_mean_swing_len = struct_walk_back_steps.left_swing_len_mean;
comb_struct.walking_left_mean_swing_duration = struct_walk_back_steps.left_swing_duration_mean;
comb_struct.walking_left_mean_swing_speed = struct_walk_back_steps.left_swing_speed_mean;
comb_struct.walking_left_mean_swing_ratio = struct_walk_back_steps.left_swing_ratio_mean;

comb_struct.walking_left_mean_stance_len = struct_walk_back_steps.left_stance_len_mean;
comb_struct.walking_left_mean_stance_duration = struct_walk_back_steps.left_stance_duration_mean;
comb_struct.walking_left_mean_stance_speed = struct_walk_back_steps.left_stance_speed_mean;
comb_struct.walking_left_mean_stance_ratio = struct_walk_back_steps.left_stance_ratio_mean;

comb_struct.walking_left_mean_temp_ratio = struct_walk_back_steps.left_temp_ratio_mean;

comb_struct.walking_left_step_count = struct_walk_back_steps.left_step_count;

comb_struct.walking_left_std_step_len = struct_walk_back_steps.left_step_len_std;
comb_struct.walking_left_std_step_duration = struct_walk_back_steps.left_step_duration_std;
comb_struct.walking_left_std_step_speed = struct_walk_back_steps.left_step_speed_std;
comb_struct.walking_left_std_step_height = struct_walk_back_steps.left_step_height_std;
comb_struct.walking_left_std_step_width = struct_walk_back_steps.left_step_width_std;

comb_struct.walking_left_std_swing_len = struct_walk_back_steps.left_swing_len_std;
comb_struct.walking_left_std_swing_duration = struct_walk_back_steps.left_swing_duration_std;
comb_struct.walking_left_std_swing_speed = struct_walk_back_steps.left_swing_speed_std;
comb_struct.walking_left_std_swing_ratio = struct_walk_back_steps.left_swing_ratio_std;

comb_struct.walking_left_std_stance_len = struct_walk_back_steps.left_stance_len_std;
comb_struct.walking_left_std_stance_duration = struct_walk_back_steps.left_stance_duration_std;
comb_struct.walking_left_std_stance_speed = struct_walk_back_steps.left_stance_speed_std;
comb_struct.walking_left_std_stance_ratio = struct_walk_back_steps.left_stance_ratio_std;

comb_struct.walking_left_std_temp_ratio = struct_walk_back_steps.left_temp_ratio_std;


comb_struct.walking_norm_left_mean_step_len = struct_walk_back_steps.norm_left_step_len_mean;
comb_struct.walking_norm_left_mean_step_duration = struct_walk_back_steps.norm_left_step_duration_mean;
comb_struct.walking_norm_left_mean_step_speed = struct_walk_back_steps.norm_left_step_speed_mean;
comb_struct.walking_norm_left_mean_swing_len = struct_walk_back_steps.norm_left_swing_len_mean;
comb_struct.walking_norm_left_mean_swing_duration = struct_walk_back_steps.norm_left_swing_duration_mean;
comb_struct.walking_norm_left_mean_swing_speed = struct_walk_back_steps.norm_left_swing_speed_mean;
comb_struct.walking_norm_left_mean_stance_len = struct_walk_back_steps.norm_left_stance_len_mean;
comb_struct.walking_norm_left_mean_stance_duration = struct_walk_back_steps.norm_left_stance_duration_mean;
comb_struct.walking_norm_left_mean_stance_speed = struct_walk_back_steps.norm_left_stance_speed_mean;

comb_struct.walking_norm_left_std_step_len = struct_walk_back_steps.norm_left_step_len_std;
comb_struct.walking_norm_left_std_step_duration = struct_walk_back_steps.norm_left_step_duration_std;
comb_struct.walking_norm_left_std_step_speed = struct_walk_back_steps.norm_left_step_speed_std;
comb_struct.walking_norm_left_std_swing_len = struct_walk_back_steps.norm_left_swing_len_std;
comb_struct.walking_norm_left_std_swing_duration = struct_walk_back_steps.norm_left_swing_duration_std;
comb_struct.walking_norm_left_std_swing_speed = struct_walk_back_steps.norm_left_swing_speed_std;
comb_struct.walking_norm_left_std_stance_len = struct_walk_back_steps.norm_left_stance_len_std;
comb_struct.walking_norm_left_std_stance_duration = struct_walk_back_steps.norm_left_stance_duration_std;
comb_struct.walking_norm_left_std_stance_speed = struct_walk_back_steps.norm_left_stance_speed_std;

clear struct_walk_back_steps
rmpath('./Walking Steps');

disp('    Walking Steps Saved');

% Walking Angles
addpath('./Walking Angles');

load(strcat(p_num,'_walk_back_angle.mat'));

comb_struct.walking_foot_slice = struct_walk_back_angle.leg;
comb_struct.walking_angle_steps = struct_walk_back_angle.steps;

comb_struct.walking_mean_transition_SHO_XY_ang = struct_walk_back_angle.transition_SHO_XY_ang_mean;
comb_struct.walking_mean_transition_SHO_XZ_ang = struct_walk_back_angle.transition_SHO_XZ_ang_mean;
comb_struct.walking_mean_transition_head_ang = struct_walk_back_angle.transition_head_ang_mean;
comb_struct.walking_mean_transition_back_ang = struct_walk_back_angle.transition_back_ang_mean;
comb_struct.walking_mean_transition_rwrst_dist = struct_walk_back_angle.transition_rwrst_dist_mean;
comb_struct.walking_mean_transition_lwrst_dist = struct_walk_back_angle.transition_lwrst_dist_mean;

comb_struct.walking_mean_std_transition_SHO_XY_ang = struct_walk_back_angle.transition_SHO_XY_ang_mean_std;
comb_struct.walking_mean_std_transition_SHO_XZ_ang = struct_walk_back_angle.transition_SHO_XZ_ang_mean_std;
comb_struct.walking_mean_std_transition_head_ang = struct_walk_back_angle.transition_head_ang_mean_std;
comb_struct.walking_mean_std_transition_back_ang = struct_walk_back_angle.transition_back_ang_mean_std;
comb_struct.walking_mean_std_transition_rwrst_dist = struct_walk_back_angle.transition_rwrst_dist_mean_std;
comb_struct.walking_mean_std_transition_lwrst_dist = struct_walk_back_angle.transition_lwrst_dist_mean_std;

comb_struct.walking_std_transition_SHO_XY_ang = struct_walk_back_angle.transition_SHO_XY_ang_std;
comb_struct.walking_std_transition_SHO_XZ_ang = struct_walk_back_angle.transition_SHO_XZ_ang_std;
comb_struct.walking_std_transition_head_ang = struct_walk_back_angle.transition_head_ang_std;
comb_struct.walking_std_transition_back_ang = struct_walk_back_angle.transition_back_ang_std;
comb_struct.walking_std_transition_rwrst_dist = struct_walk_back_angle.transition_rwrst_dist_std;
comb_struct.walking_std_transition_lwrst_dist = struct_walk_back_angle.transition_lwrst_dist_std;

comb_struct.walking_mean_SHO_XY_ang = struct_walk_back_angle.SHO_XY_ang_mean;
comb_struct.walking_mean_SHO_XZ_ang = struct_walk_back_angle.SHO_XZ_ang_mean;
comb_struct.walking_mean_head_ang = struct_walk_back_angle.head_ang_mean;
comb_struct.walking_mean_back_ang = struct_walk_back_angle.back_ang_mean;
comb_struct.walking_mean_rwrst_dist = struct_walk_back_angle.rwrst_dist_mean;
comb_struct.walking_mean_lwrst_dist = struct_walk_back_angle.lwrst_dist_mean;

comb_struct.walking_mean_std_SHO_XY_ang = struct_walk_back_angle.SHO_XY_ang_mean_std;
comb_struct.walking_mean_std_SHO_XZ_ang = struct_walk_back_angle.SHO_XZ_ang_mean_std;
comb_struct.walking_mean_std_head_ang = struct_walk_back_angle.head_ang_mean_std;
comb_struct.walking_mean_std_back_ang = struct_walk_back_angle.back_ang_mean_std;
comb_struct.walking_mean_std_rwrst_dist = struct_walk_back_angle.rwrst_dist_mean_std;
comb_struct.walking_mean_std_lwrst_dist = struct_walk_back_angle.lwrst_dist_mean_std;

comb_struct.walking_std_SHO_XY_ang = struct_walk_back_angle.SHO_XY_ang_std;
comb_struct.walking_std_SHO_XZ_ang = struct_walk_back_angle.SHO_XZ_ang_std;
comb_struct.walking_std_head_ang = struct_walk_back_angle.head_ang_std;
comb_struct.walking_std_back_ang = struct_walk_back_angle.back_ang_std;
comb_struct.walking_std_rwrst_dist = struct_walk_back_angle.rwrst_dist_std;
comb_struct.walking_std_lwrst_dist = struct_walk_back_angle.lwrst_dist_std;

clear struct_walk_back_angle
rmpath('./Walking Angles');

disp('    Walking Angles Saved');

%% SAVE .MAT FILE

if init_struct.mat_save
        
    save(strcat(p_num,'_struct.mat'),'comb_struct');
    writetable(struct2table(comb_struct), strcat(p_num,'_struct.xlsx'));
    
    disp('    .MAT + Excel Saved');

end

end
