function [struct_walk_back_steps] = walk_back_steps(init_struct,p_num,ind_end,Mrk_Data,Static_Data)
% function [struct_walk_back_steps] = walk_back_steps(init_struct,p_num,Mrk_Data)
% 
%  This function calculates the walking back steps parameters 
%  from marker data while saving struct_walk_back_steps and 
%  plots all the graph and .mat file
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Number (string)
%           ind_end - End Index to slice data (index)
%           Mrk_Data - Marker Data (table)
%           Static_Data - Static File parameters (struct)
%     
%   OUTPUT: struct_walk_back_steps - Walking back step parameters (struct)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% EXCEPTIONS

exceptions ={'1011_dy6';'1011_dy8';'1013_dy6';'1022_dy4';'1022_dy5';'1022_dy6';'1022_dy7';'1022_dy8';};

if ismember(p_num,exceptions)
    min_dist = 60;
else
    min_dist = 100;
end

t_exceptions ={'1009_dy4';'1007_dy5';'1012_dy5';'1024_dy4';};

%% CHANGE FOLDER

if ~exist('Walking Steps', 'dir') % Check if folder exist
    mkdir('Walking Steps'); % make new folder
end

cd('Walking Steps'); % Change directory to new folder

disp("Walking Steps: ");

%% GAIT INITIALIZATION DETECTION

walk_start = ind_end + 600; % ind_end = 5 sec before first step

time_frst_stp = Mrk_Data.Time(walk_start);

%% MARKER DATA

% Base of spine / Iliopelviic (IP)
mrkr_IPx = Mrk_Data.IPX(:); % IPX

% Right Ankle (RANK)
mrkr_RANKx = Mrk_Data.RANKX(:); % RANKX
mrkr_RANKy = Mrk_Data.RANKY(:); % RANKY
mrkr_RANKz = Mrk_Data.RANKZ(:); % RANKZ

% Left Ankle (LANK)
mrkr_LANKx = Mrk_Data.LANKX(:); % LANKX
mrkr_LANKy = Mrk_Data.LANKY(:); % LANKY
mrkr_LANKz = Mrk_Data.LANKZ(:); % LANKZ

% Right Heel (RHEEL)
mrkr_RHEELz = Mrk_Data.RHEELZ(:); % RHEELZ
mrkr_RHEELx = Mrk_Data.RHEELX(:); % RHEELX

% Left Heel (LHEEL)
mrkr_LHEELz = Mrk_Data.LHEELZ(:); % LHEELZ
mrkr_LHEELx = Mrk_Data.LHEELX(:); % LHEELX

% Right Toe (RTOE)
mrkr_RTOEz = Mrk_Data.RTOEZ(:); % RTOEZ
mrkr_RTOEx = Mrk_Data.RTOEX(:); % RTOEX

% Left Toe (LTOE)
mrkr_LTOEz = Mrk_Data.LTOEZ(:); % LTOEZ
mrkr_LTOEx = Mrk_Data.LTOEX(:); % LTOEX

% Filter Data
mrkr_IPx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_IPx)./1000; % m

mrkr_RANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKx)./1000; % m
mrkr_RANKy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKy)./1000; % m
mrkr_RANKz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKz)./1000; % m

mrkr_LANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKx)./1000; % m
mrkr_LANKy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKy)./1000; % m
mrkr_LANKz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKz)./1000; % m

mrkr_RHEELz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELz)./1000; % m
mrkr_RHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELx)./1000; % m

mrkr_LHEELz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELz)./1000; % m
mrkr_LHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELx)./1000; % m

mrkr_RTOEz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEz)./1000; % m
mrkr_RTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEx)./1000; % m

mrkr_LTOEz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEz)./1000; % m
mrkr_LTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEx)./1000; % m

%% WALKING BACK DETECT

% Time vector
Time_Step = Mrk_Data.Time(:);
RANKX_diff = diff(mrkr_RANKx_filter);
LANKX_diff = diff(mrkr_LANKx_filter);

if ismember(p_num,t_exceptions)
    forward_frames = 15;
    walk_start = min(find((RANKX_diff(120:end)<-0.01),1,'first'),find((LANKX_diff(120:end)<-0.01),1,'first'));
else
    forward_frames = 30;
    walk_start = max(find((RANKX_diff(120:end)<-0.01),1,'first'),find((LANKX_diff(120:end)<-0.01),1,'first'));
end

walk_start = walk_start + forward_frames;

Time_turn = Time_Step(1:end-1);
turn_time = Time_turn(walk_start);

figure;
subplot(2,1,1);
plot(Time_turn(walk_start-600:end),RANKX_diff(walk_start-600:end));
hold on;
line([turn_time turn_time],[-0.04 0.04],'Color',[1 0 0]);
title("Right Ankle");
ylabel("RANKX Diff");

subplot(2,1,2);
plot(Time_turn(walk_start-600:end),LANKX_diff(walk_start-600:end));
line([turn_time turn_time],[-0.04 0.04],'Color',[1 0 0]);
title("Left Ankle");
ylabel("LANKX Diff");

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'walk_back_detect.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Step Detction Calculated');


%% CLEAN MARKER

% Base of spine / Iliopelviic (IP)
IPX = mrkr_IPx_filter(walk_start:end); % IPX

% Right Ankle (RANK)
RANKX = mrkr_RANKx_filter(walk_start:end); % RANKX
RANKY = mrkr_RANKy_filter(walk_start:end); % RANKY
RANKZ = mrkr_RANKz_filter(walk_start:end); % RANKZ

% Left Ankle (LANK)
LANKX = mrkr_LANKx_filter(walk_start:end); % LANKX
LANKY = mrkr_LANKy_filter(walk_start:end); % LANKY
LANKZ = mrkr_LANKz_filter(walk_start:end); % LANKZ

% Right Heel (RHEEL)
RHEELZ = mrkr_RHEELz_filter(walk_start:end); % RHEELZ
RHEELX = mrkr_RHEELx_filter(walk_start:end); % RHEELX

% Left Heel (LHEEL)
LHEELZ = mrkr_LHEELz_filter(walk_start:end); % LHEELZ
LHEELX = mrkr_LHEELx_filter(walk_start:end); % LHEELX

% Right Toe (RTOE)
RTOEZ = mrkr_RTOEz_filter(walk_start:end); % RTOEZ
RTOEX = mrkr_RTOEx_filter(walk_start:end); % RTOEX

% Left Toe (LTOE)
LTOEZ = mrkr_LTOEz_filter(walk_start:end); % LTOEZ
LTOEX = mrkr_LTOEx_filter(walk_start:end); % LTOEX

clear mrkr_* % clear unnecessary variables 

%% DETECT STEPS

% Time vector
Time_Step = Mrk_Data.Time(walk_start:end);

% Distance
RHeel_stps_dist = (IPX-RHEELX);
LHeel_stps_dist = (IPX-LHEELX);

RToe_stps_dist = -(IPX-RTOEX);
LToe_stps_dist = -(IPX-LTOEX);

[~, peak_rstp_indx] = findpeaks(RHeel_stps_dist,'MinPeakDistance',min_dist);
[~, peak_lstp_indx] = findpeaks(LHeel_stps_dist,'MinPeakDistance',min_dist);

[~, valley_rstp_indx] = findpeaks(RToe_stps_dist,'MinPeakDistance',min_dist);
[~, valley_lstp_indx] = findpeaks(LToe_stps_dist,'MinPeakDistance',min_dist);

% heel strike index
rhs_indx = [peak_rstp_indx,ones(size(peak_rstp_indx))];
lhs_indx = [peak_lstp_indx,ones(size(peak_lstp_indx))*-1];
hs_indx = sortrows(vertcat(rhs_indx,lhs_indx),1);

% toe-off index
rto_indx = [valley_rstp_indx,ones(size(valley_rstp_indx))];
lto_indx = [valley_lstp_indx,ones(size(valley_lstp_indx))*-1];
to_indx = sortrows(vertcat(rto_indx,lto_indx),1);

% kepp same number of indexes
lhs_indx = hs_indx(hs_indx(:,2)==-1);
rhs_indx = hs_indx(hs_indx(:,2)==1);

lto_indx = to_indx(to_indx(:,2)==-1);
rto_indx = to_indx(to_indx(:,2)==1);

if (length(rhs_indx))~=(length(rto_indx)) && (length(rhs_indx)>1)
    min_r = min(length(rhs_indx),length(rto_indx));
    if length(rto_indx)>min_r
        len_indx = length(rto_indx)+1;
        rto_indx = rto_indx(len_indx-min_r:end);
    else
        len_indx = length(rhs_indx)+1;
        rhs_indx = rhs_indx(len_indx-min_r:end);
    end
    
elseif (length(lhs_indx))~=(length(lto_indx)) && (length(lhs_indx)>1)
    min_l = min(length(lhs_indx),length(lto_indx));
    if length(lto_indx)>min_l
        len_indx = length(lto_indx)+1;
        lto_indx = lto_indx(len_indx-min_l:end);
    else
        len_indx = length(lhs_indx)+1;
        lhs_indx = lhs_indx(len_indx-min_l:end);
    end
end

%% PLOTS

% Distance Plots

% Right Leg
figure;
subplot(2,1,1);
plot(Time_Step,RHeel_stps_dist);
hold on;
plot(Time_Step(rhs_indx),RHeel_stps_dist(rhs_indx),'*');
hold off;
title("IP to RHEEL ");

subplot(2,1,2);
plot(Time_Step,RToe_stps_dist);
hold on;
plot(Time_Step(rto_indx),RToe_stps_dist(rto_indx),'+');
hold off;
title("IP to RTOE");

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'IP_2_Right.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Left Leg
figure;
subplot(2,1,1);
plot(Time_Step,LHeel_stps_dist);
hold on;
plot(Time_Step(lhs_indx),LHeel_stps_dist(lhs_indx),'*');
hold off;
title("IP to LHEEL ");

subplot(2,1,2);
plot(Time_Step,LToe_stps_dist);
hold on;
plot(Time_Step(lto_indx),LToe_stps_dist(lto_indx),'+');
hold off;
title("IP to LTOE");

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'IP_2_Left.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Ankle Plots

% RANK 
figure;
subplot(2,1,1);
plot(Time_Step,RANKZ);
hold on;
plot(Time_Step(rhs_indx),RANKZ(rhs_indx),'g*');
plot(Time_Step(rto_indx),RANKZ(rto_indx),'r+');
hold off;
title("RANK across Recording");
ylabel('RANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(2,1,2);
plot(Time_Step,RANKX);
hold on;
plot(Time_Step(rhs_indx),RANKX(rhs_indx),'g*');
plot(Time_Step(rto_indx),RANKX(rto_indx),'r+');
hold off;
title("RANK across Recording");
ylabel('RANK X-Axis'),xlabel(strcat('Time (secs)'));

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'RANK_walk_back.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% LANK
figure;
subplot(2,1,1);
plot(Time_Step,LANKZ);
hold on;
plot(Time_Step(lhs_indx),LANKZ(lhs_indx),'g*');
plot(Time_Step(lto_indx),LANKZ(lto_indx),'r+');
hold off;
title("LANK across Recording");
ylabel('LANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(2,1,2);
plot(Time_Step,LANKX);
hold on;
plot(Time_Step(lhs_indx),LANKX(lhs_indx),'g*');
plot(Time_Step(lto_indx),LANKX(lto_indx),'r+');
hold off;
title("LANK across Recording");
ylabel('LANK X-Axis'),xlabel(strcat('Time (secs)'));


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'LANK_walk_back.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% ANKLE X-Z

figure;
subplot(2,1,1);
plot(RANKX,RANKZ);
hold on;
plot(RANKX(rhs_indx),RANKZ(rhs_indx),'g*');
plot(RANKX(rto_indx),RANKZ(rto_indx),'r+');
hold off;
title("RANK X vs Z");
ylabel('RANK X-Axis'),xlabel('RANK X-Axis');

subplot(2,1,2);
plot(LANKX,LANKZ);
hold on;
plot(LANKX(lhs_indx),LANKZ(lhs_indx),'g*');
plot(LANKX(lto_indx),LANKZ(lto_indx),'r+');
hold off;
title("LANK X vs Z");
ylabel('LANK Z-Axis'),xlabel('LANK X-Axis');


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'ANK_XZ_walk_back.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Plots Saved');

%% GAIT INITIALIZATION
struct_walk_back_steps = struct;
left_leg = struct;
right_leg = struct;

% Left Foot
j=1;
for i=2:length(lto_indx)
    left_leg(j).steps_count = length(lto_indx)-1;
    left_leg(j).step_len = norm(LANKX(lto_indx(i))-RANKX(lto_indx(i-1)));
    left_leg(j).step_duration = Time_Step(lto_indx(i))-Time_Step(lto_indx(i-1));
    left_leg(j).step_speed = left_leg(j).step_len/left_leg(j).step_duration;
    
    % Normalised Values
    left_leg(j).norm_step_len = left_leg(j).step_len/Static_Data.Leg_norm;
    left_leg(j).norm_step_duration = left_leg(j).step_duration/Static_Data.Time_norm;
    left_leg(j).norm_step_speed = left_leg(j).step_speed/Static_Data.Vel_norm;
    
    left_leg(j).step_height = max(LANKZ(lto_indx(i-1):lto_indx(i)))-min(LANKZ(lto_indx(i-1):lto_indx(i)));
    left_leg(j).step_width = norm(LANKY(lto_indx(i))-RANKY(lto_indx(i)));
    
    left_leg(j).swing_len = norm(LANKX(lhs_indx(i-1))-RANKX(lhs_indx(i-1)));
    left_leg(j).swing_duration = Time_Step(lhs_indx(i-1))-Time_Step(lto_indx(i-1));
    left_leg(j).swing_speed = left_leg(j).swing_len/left_leg(j).swing_duration;
    left_leg(j).swing_ratio = left_leg(j).swing_duration/left_leg(j).step_duration;
    
    % Normalised Values
    left_leg(j).norm_swing_len = left_leg(j).swing_len/Static_Data.Leg_norm;
    left_leg(j).norm_swing_duration = left_leg(j).swing_duration/Static_Data.Time_norm;
    left_leg(j).norm_swing_speed = left_leg(j).swing_speed/Static_Data.Vel_norm;

    left_leg(j).stance_len = norm(LANKX(lto_indx(i))-RANKX(lto_indx(i)));
    left_leg(j).stance_duration = Time_Step(lto_indx(i))-Time_Step(lhs_indx(i-1));
    left_leg(j).stance_speed = left_leg(j).stance_len/left_leg(j).stance_duration;
    left_leg(j).stance_ratio = left_leg(j).stance_duration/left_leg(j).step_duration;
    
    % Normalised Values
    left_leg(j).norm_stance_len = left_leg(j).stance_len/Static_Data.Leg_norm;
    left_leg(j).norm_stance_duration = left_leg(j).stance_duration/Static_Data.Time_norm;
    left_leg(j).norm_stance_speed = left_leg(j).stance_speed/Static_Data.Vel_norm;
    
    left_leg(j).temp_ratio = left_leg(j).stance_duration/left_leg(j).swing_duration;
    
    j=j+1;
end


% Right Foot first off
j=1;
for i=2:length(rto_indx)
    right_leg(j).steps_count = length(rto_indx)-1;
    right_leg(j).step_len = norm(RANKX(rto_indx(i))-LANKX(rto_indx(i-1)));
    right_leg(j).step_duration = Time_Step(rto_indx(i))-Time_Step(rto_indx(i-1));
    right_leg(j).step_speed = right_leg(j).step_len/right_leg(j).step_duration;
    
    % Normalised Values
    right_leg(j).norm_step_len = right_leg(j).step_len/Static_Data.Leg_norm;
    right_leg(j).norm_step_duration = right_leg(j).step_duration/Static_Data.Time_norm;
    right_leg(j).norm_step_speed = right_leg(j).step_speed/Static_Data.Vel_norm;
    
    right_leg(j).step_height = max(RANKZ(rto_indx(i-1):rto_indx(i)))-min(RANKZ(rto_indx(i-1):rto_indx(i)));
    right_leg(j).step_width = norm(RANKY(rto_indx(i))-LANKY(rto_indx(i)));
    
    right_leg(j).swing_len = norm(RANKX(rhs_indx(i-1))-LANKX(rhs_indx(i-1)));
    right_leg(j).swing_duration = Time_Step(rhs_indx(i-1))-Time_Step(rto_indx(i-1));
    right_leg(j).swing_speed = right_leg(j).swing_len/right_leg(j).swing_duration;
    right_leg(j).swing_ratio = right_leg(j).swing_duration/right_leg(j).step_duration;
    
    % Normalised Values
    right_leg(j).norm_swing_len = right_leg(j).swing_len/Static_Data.Leg_norm;
    right_leg(j).norm_swing_duration = right_leg(j).swing_duration/Static_Data.Time_norm;
    right_leg(j).norm_swing_speed = right_leg(j).swing_speed/Static_Data.Vel_norm;
    
    right_leg(j).stance_len = norm(RANKX(rto_indx(i))-LANKX(rto_indx(i)));
    right_leg(j).stance_duration = Time_Step(rto_indx(i))-Time_Step(rhs_indx(i-1));
    right_leg(j).stance_speed = right_leg(j).stance_len/right_leg(j).stance_duration;
    right_leg(j).stance_ratio = right_leg(j).stance_duration/right_leg(j).step_duration;
    
    % Normalised Values
    right_leg(j).norm_stance_len = right_leg(j).stance_len/Static_Data.Leg_norm;
    right_leg(j).norm_stance_duration = right_leg(j).stance_duration/Static_Data.Time_norm;
    right_leg(j).norm_stance_speed = right_leg(j).stance_speed/Static_Data.Vel_norm;
    
    right_leg(j).temp_ratio = right_leg(j).stance_duration/right_leg(j).swing_duration;
    
    j=j+1;
end

struct_walk_back_steps.right_leg = right_leg;
struct_walk_back_steps.left_leg = left_leg;

struct_walk_back_steps.right_step_count = right_leg(1).steps_count;
% Right Leg Mean
struct_walk_back_steps.right_step_len_mean = mean([right_leg.step_len]);
struct_walk_back_steps.right_step_duration_mean = mean([right_leg.step_duration]);
struct_walk_back_steps.right_step_speed_mean = mean([right_leg.step_speed]);
struct_walk_back_steps.right_step_height_mean = mean([right_leg.step_height]);
struct_walk_back_steps.right_step_width_mean = mean([right_leg.step_width]);
struct_walk_back_steps.right_swing_len_mean = mean([right_leg.swing_len]);
struct_walk_back_steps.right_swing_duration_mean = mean([right_leg.swing_duration]);
struct_walk_back_steps.right_swing_speed_mean = mean([right_leg.swing_speed]);
struct_walk_back_steps.right_swing_ratio_mean = mean([right_leg.swing_ratio]);
struct_walk_back_steps.right_stance_len_mean = mean([right_leg.stance_len]);
struct_walk_back_steps.right_stance_duration_mean = mean([right_leg.stance_duration]);
struct_walk_back_steps.right_stance_speed_mean = mean([right_leg.stance_speed]);
struct_walk_back_steps.right_stance_ratio_mean = mean([right_leg.stance_ratio]);
struct_walk_back_steps.right_temp_ratio_mean = mean([right_leg.temp_ratio]);
% Right Leg Normalized Mean
struct_walk_back_steps.norm_right_step_len_mean = mean([right_leg.norm_step_len]);
struct_walk_back_steps.norm_right_step_duration_mean = mean([right_leg.norm_step_duration]);
struct_walk_back_steps.norm_right_step_speed_mean = mean([right_leg.norm_step_speed]);
struct_walk_back_steps.norm_right_swing_len_mean = mean([right_leg.norm_swing_len]);
struct_walk_back_steps.norm_right_swing_duration_mean = mean([right_leg.norm_swing_duration]);
struct_walk_back_steps.norm_right_swing_speed_mean = mean([right_leg.norm_swing_speed]);
struct_walk_back_steps.norm_right_stance_len_mean = mean([right_leg.norm_stance_len]);
struct_walk_back_steps.norm_right_stance_duration_mean = mean([right_leg.norm_stance_duration]);
struct_walk_back_steps.norm_right_stance_speed_mean = mean([right_leg.norm_stance_speed]);
% Right Leg Std
struct_walk_back_steps.right_step_len_std = std([right_leg.step_len]);
struct_walk_back_steps.right_step_duration_std = std([right_leg.step_duration]);
struct_walk_back_steps.right_step_speed_std = std([right_leg.step_speed]);
struct_walk_back_steps.right_step_height_std = std([right_leg.step_height]);
struct_walk_back_steps.right_step_width_std = std([right_leg.step_width]);
struct_walk_back_steps.right_swing_len_std = std([right_leg.swing_len]);
struct_walk_back_steps.right_swing_duration_std = std([right_leg.swing_duration]);
struct_walk_back_steps.right_swing_speed_std = std([right_leg.swing_speed]);
struct_walk_back_steps.right_swing_ratio_std = std([right_leg.swing_ratio]);
struct_walk_back_steps.right_stance_len_std = std([right_leg.stance_len]);
struct_walk_back_steps.right_stance_duration_std = std([right_leg.stance_duration]);
struct_walk_back_steps.right_stance_speed_std = std([right_leg.stance_speed]);
struct_walk_back_steps.right_stance_ratio_std = std([right_leg.stance_ratio]);
struct_walk_back_steps.right_temp_ratio_std = std([right_leg.temp_ratio]);
% Right Leg Normalized Std
struct_walk_back_steps.norm_right_step_len_std = std([right_leg.norm_step_len]);
struct_walk_back_steps.norm_right_step_duration_std = std([right_leg.norm_step_duration]);
struct_walk_back_steps.norm_right_step_speed_std = std([right_leg.norm_step_speed]);
struct_walk_back_steps.norm_right_swing_len_std = std([right_leg.norm_swing_len]);
struct_walk_back_steps.norm_right_swing_duration_std = std([right_leg.norm_swing_duration]);
struct_walk_back_steps.norm_right_swing_speed_std = std([right_leg.norm_swing_speed]);
struct_walk_back_steps.norm_right_stance_len_std = std([right_leg.norm_stance_len]);
struct_walk_back_steps.norm_right_stance_duration_std = std([right_leg.norm_stance_duration]);
struct_walk_back_steps.norm_right_stance_speed_std = std([right_leg.norm_stance_speed]);

struct_walk_back_steps.left_step_count = left_leg(1).steps_count;
% Left Leg Mean
struct_walk_back_steps.left_step_len_mean = mean([left_leg.step_len]);
struct_walk_back_steps.left_step_duration_mean = mean([left_leg.step_duration]);
struct_walk_back_steps.left_step_speed_mean = mean([left_leg.step_speed]);
struct_walk_back_steps.left_step_height_mean = mean([left_leg.step_height]);
struct_walk_back_steps.left_step_width_mean = mean([left_leg.step_width]);
struct_walk_back_steps.left_swing_len_mean = mean([left_leg.swing_len]);
struct_walk_back_steps.left_swing_duration_mean = mean([left_leg.swing_duration]);
struct_walk_back_steps.left_swing_speed_mean = mean([left_leg.swing_speed]);
struct_walk_back_steps.left_swing_ratio_mean = mean([left_leg.swing_ratio]);
struct_walk_back_steps.left_stance_len_mean = mean([left_leg.stance_len]);
struct_walk_back_steps.left_stance_duration_mean = mean([left_leg.stance_duration]);
struct_walk_back_steps.left_stance_speed_mean = mean([left_leg.stance_speed]);
struct_walk_back_steps.left_stance_ratio_mean = mean([left_leg.stance_ratio]);
struct_walk_back_steps.left_temp_ratio_mean = mean([left_leg.temp_ratio]);
% Left Leg Normalized Mean
struct_walk_back_steps.norm_left_step_len_mean = mean([left_leg.norm_step_len]);
struct_walk_back_steps.norm_left_step_duration_mean = mean([left_leg.norm_step_duration]);
struct_walk_back_steps.norm_left_step_speed_mean = mean([left_leg.norm_step_speed]);
struct_walk_back_steps.norm_left_swing_len_mean = mean([left_leg.norm_swing_len]);
struct_walk_back_steps.norm_left_swing_duration_mean = mean([left_leg.norm_swing_duration]);
struct_walk_back_steps.norm_left_swing_speed_mean = mean([left_leg.norm_swing_speed]);
struct_walk_back_steps.norm_left_stance_len_mean = mean([left_leg.norm_stance_len]);
struct_walk_back_steps.norm_left_stance_duration_mean = mean([left_leg.norm_stance_duration]);
struct_walk_back_steps.norm_left_stance_speed_mean = mean([left_leg.norm_stance_speed]);
% Left Leg Std
struct_walk_back_steps.left_step_len_std = std([left_leg.step_len]);
struct_walk_back_steps.left_step_duration_std = std([left_leg.step_duration]);
struct_walk_back_steps.left_step_speed_std = std([left_leg.step_speed]);
struct_walk_back_steps.left_step_height_std = std([left_leg.step_height]);
struct_walk_back_steps.left_step_width_std = std([left_leg.step_width]);
struct_walk_back_steps.left_swing_len_std = std([left_leg.swing_len]);
struct_walk_back_steps.left_swing_duration_std = std([left_leg.swing_duration]);
struct_walk_back_steps.left_swing_speed_std = std([left_leg.swing_speed]);
struct_walk_back_steps.left_swing_ratio_std = std([left_leg.swing_ratio]);
struct_walk_back_steps.left_stance_len_std = std([left_leg.stance_len]);
struct_walk_back_steps.left_stance_duration_std = std([left_leg.stance_duration]);
struct_walk_back_steps.left_stance_speed_std = std([left_leg.stance_speed]);
struct_walk_back_steps.left_stance_ratio_std = std([left_leg.stance_ratio]);
struct_walk_back_steps.left_temp_ratio_std = std([left_leg.temp_ratio]);
% Left Leg Normalized Std
struct_walk_back_steps.norm_left_step_len_std = std([left_leg.norm_step_len]);
struct_walk_back_steps.norm_left_step_duration_std = std([left_leg.norm_step_duration]);
struct_walk_back_steps.norm_left_step_speed_std = std([left_leg.norm_step_speed]);
struct_walk_back_steps.norm_left_swing_len_std = std([left_leg.norm_swing_len]);
struct_walk_back_steps.norm_left_swing_duration_std = std([left_leg.norm_swing_duration]);
struct_walk_back_steps.norm_left_swing_speed_std = std([left_leg.norm_swing_speed]);
struct_walk_back_steps.norm_left_stance_len_std = std([left_leg.norm_stance_len]);
struct_walk_back_steps.norm_left_stance_duration_std = std([left_leg.norm_stance_duration]);
struct_walk_back_steps.norm_left_stance_speed_std = std([left_leg.norm_stance_speed]);

disp('    Step Parameters Calculated');

%% SAVE .MAT FILE

if init_struct.mat_save
    
    save(strcat(p_num,'_walk_back_steps.mat'),'struct_walk_back_steps');
    disp('    .MAT Saved');

end

%% CHANGE BACK TO DIRECTORY

cd('..')

end