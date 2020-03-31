function [struct_gait_init] = gait_init_calc(init_struct,p_num,ind_end,FP_data,Mrk_Data,Static_Data)
% function [struct_gait_init] = gait_init_detect(init_struct,p_num,FP_data,Mrk_Data)
% 
%  This function calculates the gait initialization (first step) from 
%  forceplate and marker data and calculates the gait initizialization 
%  parameters while saving struct_gait_init and plots all the graph and .mat file
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Number (string)
%           ind_end - End Index to slice data (index)
%           FP_data - Force plate data (table)
%           Mrk_Data - Marker Data (table)
%           Static_Data - Static File parameters (struct)
%     
%   OUTPUT: struct_gait_init - Gait Initialization parameters (struct)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% CHANGE FOLDER

if ~exist('Gait Init', 'dir') % Check if folder exist
    mkdir('Gait Init'); % make new folder
end

cd('Gait Init'); % Change directory to new folder

disp("Gait Initilization: ");

%% GAIT INITIALIZATION DETECTION

indx_frst_stp = ind_end + 600; % ind_end = 5 sec before first step

gait_init_start = indx_frst_stp - 180;
gait_init_end = indx_frst_stp + 480; % 4 sec after step

time_frst_stp = Mrk_Data.Time(indx_frst_stp);

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

%% CLEAN MARKER DATA

% Base of spine / Iliopelviic (IP)
IPX = mrkr_IPx_filter(gait_init_start:gait_init_end); % IPX

% Right Ankle (RANK)
RANKX = mrkr_RANKx_filter(gait_init_start:gait_init_end); % RANKX
RANKY = mrkr_RANKy_filter(gait_init_start:gait_init_end); % RANKY
RANKZ = mrkr_RANKz_filter(gait_init_start:gait_init_end); % RANKZ

% Left Ankle (LANK)
LANKX = mrkr_LANKx_filter(gait_init_start:gait_init_end); % LANKX
LANKY = mrkr_LANKy_filter(gait_init_start:gait_init_end); % LANKY
LANKZ = mrkr_LANKz_filter(gait_init_start:gait_init_end); % LANKZ

% Right Heel (RHEEL)
RHEELZ = mrkr_RHEELz_filter(gait_init_start:gait_init_end); % RHEELZ
RHEELX = mrkr_RHEELx_filter(gait_init_start:gait_init_end); % RHEELX

% Left Heel (LHEEL)
LHEELZ = mrkr_LHEELz_filter(gait_init_start:gait_init_end); % LHEELZ
LHEELX = mrkr_LHEELx_filter(gait_init_start:gait_init_end); % LHEELX

% Right Toe (RTOE)
RTOEZ = mrkr_RTOEz_filter(gait_init_start:gait_init_end); % RTOEZ
RTOEX = mrkr_RTOEx_filter(gait_init_start:gait_init_end); % RTOEX

% Left Toe (LTOE)
LTOEZ = mrkr_LTOEz_filter(gait_init_start:gait_init_end); % LTOEZ
LTOEX = mrkr_LTOEx_filter(gait_init_start:gait_init_end); % LTOEX

clear mrkg*

%% DETECT STEPS

% Time vector
Time_Step = Mrk_Data.Time(gait_init_start:gait_init_end);
off_plate_indx = find(Time_Step==time_frst_stp,1,'first');

% Distance
RHeel_stps_dist = -(IPX-RHEELX);
LHeel_stps_dist = -(IPX-LHEELX);

RToe_stps_dist = (IPX-RTOEX);
LToe_stps_dist = (IPX-LTOEX);

[~, peak_rstp_indx] = findpeaks(RHeel_stps_dist,'MinPeakDistance',120);
[~, peak_lstp_indx] = findpeaks(LHeel_stps_dist,'MinPeakDistance',120);

[~, valley_rstp_indx] = findpeaks(RToe_stps_dist,'MinPeakDistance',120);
[~, valley_lstp_indx] = findpeaks(LToe_stps_dist,'MinPeakDistance',120);


% heel strike index
rhs_indx = [peak_rstp_indx,ones(size(peak_rstp_indx))];
lhs_indx = [peak_lstp_indx,ones(size(peak_lstp_indx))*-1];
hs_indx = sortrows(vertcat(rhs_indx,lhs_indx),1);

% toe-off index
rto_indx = [valley_rstp_indx,ones(size(valley_rstp_indx))];
lto_indx = [valley_lstp_indx,ones(size(valley_lstp_indx))*-1];
to_indx = sortrows(vertcat(rto_indx,lto_indx),1);

% kepp only first step index
[~,diff_indx]=min(abs(to_indx(:,1)-off_plate_indx));
foot_off = to_indx(diff_indx,:);
diff_indx = diff_indx+2;


to_indx = to_indx(1:diff_indx,:);
lto_indx = to_indx(to_indx(:,2)==-1);
rto_indx = to_indx(to_indx(:,2)==1);

if ismember(foot_off(:,1),rto_indx(:,1))
    rto_indx = foot_off(:,1);
    lto_indx = lto_indx(lto_indx(:,1)<rto_indx(:,1));
else
    lto_indx = foot_off(:,1);
    rto_indx = rto_indx(rto_indx(:,1)<lto_indx(:,1));
end

lto_indx = lto_indx(1);
rto_indx = rto_indx(1);

hs_indx = hs_indx(1:diff_indx,:);
lhs_indx = hs_indx(hs_indx(:,2)==-1 & hs_indx(:,1)>lto_indx);
lhs_indx = lhs_indx(1);
rhs_indx = hs_indx(hs_indx(:,2)==1 & hs_indx(:,1)>rto_indx);
rhs_indx = rhs_indx(1);

%% PLOTS

% Distance Plots

% Right Leg
figure;
subplot(2,1,1);
plot(Time_Step,RHeel_stps_dist);
hold on;
line([time_frst_stp time_frst_stp],[-0.5 0.5],'Color',[1 0 0]);
plot(Time_Step(rhs_indx),RHeel_stps_dist(rhs_indx),'*');
hold off;
title("IP to RHEEL ");

subplot(2,1,2);
plot(Time_Step,RToe_stps_dist);
hold on;
line([time_frst_stp time_frst_stp],[-1 1],'Color',[1 0 0]);
plot(Time_Step(rto_indx),RToe_stps_dist(rto_indx),'+');
hold off;
title("IP to RTOE");

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_IP_2_Right.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Left Leg
figure;
subplot(2,1,1);
plot(Time_Step,LHeel_stps_dist);
hold on;
line([time_frst_stp time_frst_stp],[-0.5 0.5],'Color',[1 0 0]);
plot(Time_Step(lhs_indx),LHeel_stps_dist(lhs_indx),'*');
hold off;
title("IP to LHEEL ");

subplot(2,1,2);
plot(Time_Step,LToe_stps_dist);
hold on;
line([time_frst_stp time_frst_stp],[-1 1],'Color',[1 0 0]);
plot(Time_Step(lto_indx),LToe_stps_dist(lto_indx),'+');
hold off;
title("IP to LTOE");

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_IP_2_Left.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Ankle Plots

% RANK 
figure;
subplot(3,1,1);
plot(FP_data.TIME(gait_init_start:gait_init_end),FP_data.Force_Z(gait_init_start:gait_init_end));
line([time_frst_stp time_frst_stp],[-1000 1000],'Color',[1 0 0])
title("Force Z across Recording")
ylabel('Force Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,2);
plot(Time_Step,RANKZ);
hold on;
line([time_frst_stp time_frst_stp],[0 0.5],'Color',[1 0 0]);
plot(Time_Step(rhs_indx),RANKZ(rhs_indx),'g*');
plot(Time_Step(rto_indx),RANKZ(rto_indx),'r+');
hold off;
title("RANK across Recording");
ylabel('RANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,3);
plot(Time_Step,RANKX);
hold on;
line([time_frst_stp time_frst_stp],[-5 5],'Color',[1 0 0]);
plot(Time_Step(rhs_indx),RANKX(rhs_indx),'g*');
plot(Time_Step(rto_indx),RANKX(rto_indx),'r+');
hold off;
title("RANK across Recording");
ylabel('RANK X-Axis'),xlabel(strcat('Time (secs)'));

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_RANK_gait_init.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% LANK
figure;
subplot(3,1,1);
plot(FP_data.TIME(gait_init_start:gait_init_end),FP_data.Force_Z(gait_init_start:gait_init_end));
line([time_frst_stp time_frst_stp],[-1000 1000],'Color',[1 0 0])
title("Force Z across Recording")
ylabel('Force Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,2);
plot(Time_Step,LANKZ);
hold on;
line([time_frst_stp time_frst_stp],[0 0.5],'Color',[1 0 0]);
plot(Time_Step(lhs_indx),LANKZ(lhs_indx),'g*');
plot(Time_Step(lto_indx),LANKZ(lto_indx),'r+');
hold off;
title("LANK across Recording");
ylabel('LANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,3);
plot(Time_Step,LANKX);
hold on;
line([time_frst_stp time_frst_stp],[-5 5],'Color',[1 0 0]);
plot(Time_Step(lhs_indx),LANKX(lhs_indx),'g*');
plot(Time_Step(lto_indx),LANKX(lto_indx),'r+');
hold off;
title("LANK across Recording");
ylabel('LANK X-Axis'),xlabel(strcat('Time (secs)'));


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_LANK_gait_init.jpg'));
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
    saveas(gcf,strcat(p_num,'_ANK_XZ_gait_init.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Plots Saved');

%% GAIT INITIALIZATION
struct_gait_init = struct;

if foot_off(:,2) == 1 % Left Foot first off
    struct_gait_init.foot_gait_init = -1;
      
    struct_gait_init.swing_len_gait_init = norm(LANKX(lhs_indx)-RANKX(lto_indx));
    struct_gait_init.swing_duration_gait_init = Time_Step(lhs_indx)-Time_Step(lto_indx);
    struct_gait_init.swing_speed_gait_init = struct_gait_init.swing_len_gait_init/struct_gait_init.swing_duration_gait_init;
    
    struct_gait_init.swing_height_gait_init = max(LANKZ(lto_indx:lhs_indx))-min(LANKZ(lto_indx:lhs_indx));
    struct_gait_init.swing_width_gait_init = norm(LANKY(lhs_indx)-RANKY(lhs_indx));
    
else % Right Foot first off
    struct_gait_init.foot_gait_init = 1;  
    
    struct_gait_init.swing_len_gait_init = norm(RANKX(rhs_indx)-LANKX(rto_indx));
    struct_gait_init.swing_duration_gait_init = Time_Step(rhs_indx)-Time_Step(rto_indx);
    struct_gait_init.swing_speed_gait_init = struct_gait_init.swing_len_gait_init/struct_gait_init.swing_duration_gait_init;
    
    struct_gait_init.swing_height_gait_init = max(RANKZ(rto_indx:rhs_indx))-min(RANKZ(rto_indx:rhs_indx));
    struct_gait_init.swing_width_gait_init = norm(RANKY(rhs_indx)-LANKY(rhs_indx));
    
end

% Normalised Values
struct_gait_init.norm_swing_len_gait_init = struct_gait_init.swing_len_gait_init/Static_Data.Leg_norm;
struct_gait_init.norm_swing_duration_gait_init = struct_gait_init.swing_duration_gait_init/Static_Data.Time_norm;
struct_gait_init.norm_swing_speed_gait_init = struct_gait_init.swing_speed_gait_init/Static_Data.Vel_norm;

disp('    Gait Initialization Calculated');

%% FIRST STEP PARAMETERS

if foot_off(:,2) == 1 % Right foot off FP
    struct_gait_init.first_step = 1;
    
    struct_gait_init.swing_len_first_step = norm(RANKX(rhs_indx)-LANKX(rto_indx));
    struct_gait_init.swing_duration_first_step = Time_Step(rhs_indx)-Time_Step(rto_indx);
    struct_gait_init.swing_speed_first_step = struct_gait_init.swing_len_first_step/struct_gait_init.swing_duration_first_step;
    
    struct_gait_init.swing_height_first_step = max(RANKZ(rto_indx:rhs_indx))-min(RANKZ(rto_indx:rhs_indx));
    struct_gait_init.swing_width_first_step = norm(RANKY(rhs_indx)-LANKY(rhs_indx));
    
else % Left foot off FP
    struct_gait_init.first_step = -1;
    
    struct_gait_init.swing_len_first_step = norm(LANKX(lhs_indx)-RANKX(lto_indx));
    struct_gait_init.swing_duration_first_step = Time_Step(lhs_indx)-Time_Step(lto_indx);
    struct_gait_init.swing_speed_first_step = struct_gait_init.swing_len_first_step/struct_gait_init.swing_duration_first_step;
    
    struct_gait_init.swing_height_first_step = max(RANKZ(lto_indx:lhs_indx))-min(RANKZ(lto_indx:lhs_indx));
    struct_gait_init.swing_width_first_step = norm(RANKY(lhs_indx)-LANKY(lhs_indx));
    
end

% Normalised Values
struct_gait_init.norm_swing_len_first_step = struct_gait_init.swing_len_first_step/Static_Data.Leg_norm;
struct_gait_init.norm_swing_duration_first_step = struct_gait_init.swing_duration_first_step/Static_Data.Time_norm;
struct_gait_init.norm_swing_speed_first_step = struct_gait_init.swing_speed_first_step/Static_Data.Vel_norm;

disp('    First Step Calculated');

%% MAT FILE

if init_struct.mat_save

    save(strcat(p_num,'_gait_init.mat'),'struct_gait_init');
    disp("    .MAT File Saved");

end
%% CHANGE BACK TO DIRECTORY

cd('..')

end