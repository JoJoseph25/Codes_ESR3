function [struct_walk_back_angle] = walk_back_angle(init_struct,p_num,ind_end,Mrk_Data)
% function [struct_walk_back_angle] = walk_back_angle(init_struct,p_num,Mrk_Data)
% 
%  This function calculates the walking back angle parameters 
%  from marker data while saving struct_walk_back_steps and 
%  plots all the graph and .mat file
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Number (string)
%           ind_end - End Index to slice data (index)
%           Mrk_Data - Marker Data (table)
%     
%   OUTPUT: struct_walk_back_angle - Walking back angle parameters (struct)
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

if ~exist('Walking Angles', 'dir') % Check if folder exist
    mkdir('Walking Angles'); % make new folder
end

cd('Walking Angles'); % Change directory to new folder

disp("Walking Angles: ");


%% GAIT INITIALIZATION DETECTION

walk_start = ind_end + 600; % ind_end = 5 sec before first step

time_frst_stp = Mrk_Data.Time(walk_start);

%% MARKER DATA

% Top of Head (TPHD)
mrkr_TPHDx = Mrk_Data.TPHDX(:); % TPHDX
mrkr_TPHDy = Mrk_Data.TPHDY(:); % TPHDY
mrkr_TPHDz = Mrk_Data.TPHDZ(:); % TPHDZ

% Base of Neck / Cervical 7 (C7)
mrkr_C7x = Mrk_Data.C7X(:); % C7X
mrkr_C7y = Mrk_Data.C7Y(:); % C7Y
mrkr_C7z = Mrk_Data.C7Z(:); % C7Z

% Right Shoulder (RSHO)
mrkr_RSHOx = Mrk_Data.RSHOX(:); % RSHOX
mrkr_RSHOy = Mrk_Data.RSHOY(:); % RSHOY
mrkr_RSHOz = Mrk_Data.RSHOZ(:); % RSHOZ

% Left Shoulder (SHO)
mrkr_LSHOx = Mrk_Data.LSHOX(:); % LSHOX
mrkr_LSHOy = Mrk_Data.LSHOY(:); % LSHOY
mrkr_LSHOz = Mrk_Data.LSHOZ(:); % LSHOZ

% Base of spine / Iliopelviic (IP)
mrkr_IPx = Mrk_Data.IPX(:); % IPX
mrkr_IPy = Mrk_Data.IPY(:); % IPY
mrkr_IPz = Mrk_Data.IPZ(:); % IPZ

% Right Wrist (RWRST)
mrkr_RWRSTz = Mrk_Data.RWRSTZ(:); % RWRSTZ
mrkr_RWRSTy = Mrk_Data.RWRSTY(:); % RWRSTX

% Left Wrist (LWRST)
mrkr_LWRSTz = Mrk_Data.LWRSTZ(:); % LWRSTZ
mrkr_LWRSTy = Mrk_Data.LWRSTY(:); % LWRSTX

% Left HIP (LHIP)
mrkr_LHIPy = Mrk_Data.LHIPY(:); % LHIPY

% Right HIP (LHIP)
mrkr_RHIPy = Mrk_Data.RHIPY(:); % RHIPY

% Right Ankle (RANK)
mrkr_RANKx = Mrk_Data.RANKX(:); % RANKX
mrkr_RANKy = Mrk_Data.RANKY(:); % RANKY
mrkr_RANKz = Mrk_Data.RANKZ(:); % RANKZ

% Left Ankle (LANK)
mrkr_LANKx = Mrk_Data.LANKX(:); % LANKX
mrkr_LANKy = Mrk_Data.LANKY(:); % LANKY
mrkr_LANKz = Mrk_Data.LANKZ(:); % LANKZ

% Right Heel (RHEEL)
mrkr_RHEELx = Mrk_Data.RHEELX(:); % RHEELX
mrkr_RHEELy = Mrk_Data.RHEELY(:); % RHEELY
mrkr_RHEELz = Mrk_Data.RHEELZ(:); % RHEELZ

% Left Heel (LHEEL)
mrkr_LHEELx = Mrk_Data.LHEELX(:); % LHEELX
mrkr_LHEELy = Mrk_Data.LHEELY(:); % LHEELY
mrkr_LHEELz = Mrk_Data.LHEELZ(:); % LHEELZ

% Right Toe (RTOE)
mrkr_RTOEx = Mrk_Data.RTOEX(:); % RTOEX
mrkr_RTOEy = Mrk_Data.RTOEY(:); % RTOEY
mrkr_RTOEz = Mrk_Data.RTOEZ(:); % RTOEZ

% Left Toe (LHEEL)
mrkr_LTOEx = Mrk_Data.LTOEX(:); % LTOEX
mrkr_LTOEy = Mrk_Data.LTOEY(:); % LTOEY
mrkr_LTOEz = Mrk_Data.LTOEZ(:); % LTOEZ

% Filter Data
mrkr_TPHDx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_TPHDx)/1000; % m
mrkr_TPHDy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_TPHDy)/1000; % m 
mrkr_TPHDz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_TPHDz)/1000; % m

mrkr_C7x_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_C7x)/1000; % m
mrkr_C7y_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_C7y)/1000; % m
mrkr_C7z_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_C7z)/1000; % m

mrkr_RSHOx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RSHOx)/1000; % m
mrkr_RSHOy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RSHOy)/1000; % m
mrkr_RSHOz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RSHOz)/1000; % m

mrkr_LSHOx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LSHOx)/1000; % m
mrkr_LSHOy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LSHOy)/1000; % m
mrkr_LSHOz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LSHOz)/1000; % m

mrkr_IPx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_IPx)/1000; % m
mrkr_IPy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_IPy)/1000; % m
mrkr_IPz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_IPz)/1000; % m

mrkr_RWRSTz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RWRSTz)/1000; % m
mrkr_LWRSTz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LWRSTz)/1000; % m

mrkr_RWRSTy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RWRSTy)/1000; % m
mrkr_LWRSTy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LWRSTy)/1000; % m

mrkr_LHIPy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHIPy)/1000; % m
mrkr_RHIPy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHIPy)/1000; % m

mrkr_RANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKx)./1000; % m
mrkr_RANKy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKy)./1000; % m
mrkr_RANKz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKz)./1000; % m

mrkr_LANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKx)./1000; % m
mrkr_LANKy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKy)./1000; % m
mrkr_LANKz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKz)./1000; % m

mrkr_RHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELx)./1000; % m
mrkr_RHEELy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELy)./1000; % m
mrkr_RHEELz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELz)./1000; % m

mrkr_LHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELx)./1000; % m
mrkr_LHEELy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELy)./1000; % m
mrkr_LHEELz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELz)./1000; % m

mrkr_RTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEx)./1000; % m
mrkr_RTOEy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEy)./1000; % m
mrkr_RTOEz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEz)./1000; % m

mrkr_LTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEx)./1000; % m
mrkr_LTOEy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEy)./1000; % m
mrkr_LTOEz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEz)./1000; % m

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

walk_start = max(find((RANKX_diff(120:end)<-0.01),1,'first'),find((LANKX_diff(120:end)<-0.01),1,'first'));
% walk_start = walk_start - 60;
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

% Top of Head (TPHD)
walk_TPHDx = mrkr_TPHDx_filter(walk_start:end);
walk_TPHDy = mrkr_TPHDy_filter(walk_start:end);
walk_TPHDz = mrkr_TPHDz_filter(walk_start:end);

% Base of Neck / Cervical 7 (C7)
walk_C7x = mrkr_C7x_filter(walk_start:end);
walk_C7y = mrkr_C7y_filter(walk_start:end);
walk_C7z = mrkr_C7z_filter(walk_start:end);

% Right Shoulder (RSHO)
walk_RSHOx = mrkr_RSHOx_filter(walk_start:end);
walk_RSHOy = mrkr_RSHOy_filter(walk_start:end);
walk_RSHOz = mrkr_RSHOz_filter(walk_start:end);

% Left Shoulder (SHO)
walk_LSHOx = mrkr_LSHOx_filter(walk_start:end);
walk_LSHOy = mrkr_LSHOy_filter(walk_start:end);
walk_LSHOz = mrkr_LSHOz_filter(walk_start:end);

% Base of spine / Iliopelviic (IP)
walk_IPx = mrkr_IPx_filter(walk_start:end);
walk_IPy = mrkr_IPy_filter(walk_start:end);
walk_IPz = mrkr_IPz_filter(walk_start:end);

% Right Wrist (RWRST)
walk_RWRSTz = mrkr_RWRSTz_filter(walk_start:end);
walk_RWRSTy = mrkr_RWRSTy_filter(walk_start:end);

% Left Wrist (LWRST)
walk_LWRSTz = mrkr_LWRSTz_filter(walk_start:end);
walk_LWRSTy = mrkr_LWRSTy_filter(walk_start:end);

% Left HIP (LHIP)
walk_LHIPy = mrkr_LHIPy_filter(walk_start:end);

% Right HIP (LHIP)
walk_RHIPy = mrkr_RHIPy_filter(walk_start:end);

% Right Ankle (RANK)
walk_RANKx = mrkr_RANKx_filter(walk_start:end);
walk_RANKy = mrkr_RANKy_filter(walk_start:end);
walk_RANKz = mrkr_RANKz_filter(walk_start:end);

% Left Ankle (LANK)
walk_LANKx = mrkr_LANKx_filter(walk_start:end);
walk_LANKy = mrkr_LANKy_filter(walk_start:end);
walk_LANKz = mrkr_LANKz_filter(walk_start:end);

% Right Heel (RHEEL)
walk_RHEELx = mrkr_RHEELx_filter(walk_start:end);
walk_RHEELy = mrkr_RHEELy_filter(walk_start:end);
walk_RHEELz = mrkr_RHEELz_filter(walk_start:end);

% Left Heel (LHEEL)
walk_LHEELx = mrkr_LHEELx_filter(walk_start:end);
walk_LHEELy = mrkr_LHEELy_filter(walk_start:end);
walk_LHEELz = mrkr_LHEELz_filter(walk_start:end);

% Right Toe (RTOE)
walk_RTOEx = mrkr_RTOEx_filter(walk_start:end);
walk_RTOEy = mrkr_RTOEy_filter(walk_start:end);
walk_RTOEz = mrkr_RTOEz_filter(walk_start:end);

% Left Toe (LTOE)
walk_LTOEx = mrkr_LTOEx_filter(walk_start:end);
walk_LTOEy = mrkr_LTOEy_filter(walk_start:end);
walk_LTOEz = mrkr_LTOEz_filter(walk_start:end);


clear mrkr_* % clear unnecessary variables 

%% DETECT STEPS

% Time vector
Time_Step = Mrk_Data.Time(walk_start:end);

% Distance
RHeel_stps_dist = (walk_IPx-walk_RHEELx);
LHeel_stps_dist = (walk_IPx-walk_LHEELx);

RToe_stps_dist = -(walk_IPx-walk_RTOEx);
LToe_stps_dist = -(walk_IPx-walk_LTOEx);

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

% Ankle Plots

% RANK 
figure;
subplot(2,1,1);
plot(Time_Step,walk_RANKz);
hold on;
plot(Time_Step(rhs_indx),walk_RANKz(rhs_indx),'g*');
plot(Time_Step(rto_indx),walk_RANKz(rto_indx),'r+');
hold off;
title("RANK across Recording");
ylabel('RANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(2,1,2);
plot(Time_Step,walk_RANKx);
hold on;
plot(Time_Step(rhs_indx),walk_RANKx(rhs_indx),'g*');
plot(Time_Step(rto_indx),walk_RANKx(rto_indx),'r+');
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
plot(Time_Step,walk_LANKz);
hold on;
plot(Time_Step(lhs_indx),walk_LANKz(lhs_indx),'g*');
plot(Time_Step(lto_indx),walk_LANKz(lto_indx),'r+');
hold off;
title("LANK across Recording");
ylabel('LANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(2,1,2);
plot(Time_Step,walk_LANKx);
hold on;
plot(Time_Step(lhs_indx),walk_LANKx(lhs_indx),'g*');
plot(Time_Step(lto_indx),walk_LANKx(lto_indx),'r+');
hold off;
title("LANK across Recording");
ylabel('LANK X-Axis'),xlabel(strcat('Time (secs)'));


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'LANK_walk_back.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Plots Saved');

%% Angle Calculation
walk_RSHO_xy = [walk_RSHOx-walk_C7x walk_RSHOy-walk_C7y];
walk_LSHO_xy = [walk_LSHOx-walk_C7x walk_LSHOy-walk_C7y];

walk_RSHO_xz = [walk_RSHOx-walk_C7x walk_RSHOz-walk_C7z];
walk_LSHO_xz = [walk_LSHOx-walk_C7x walk_LSHOz-walk_C7z];

walk_head_v1 = [walk_TPHDx-walk_C7x walk_TPHDz-walk_C7z];
walk_head_v2 = [(walk_TPHDx-walk_C7x)*0 walk_TPHDz-walk_C7z];

walk_back_v1 = [walk_IPx-walk_C7x walk_IPz-walk_C7z];
walk_back_v2 = [(walk_IPx-walk_C7x)*0 walk_IPz-walk_C7z];

% Upper Body Angles
walk_LSHO_xy_ang = atan2( walk_LSHO_xy(:,2), walk_LSHO_xy(:,1)).*180./pi; % LSHO C7 ANG
walk_RSHO_xy_ang = atan2( walk_RSHO_xy(:,2), walk_RSHO_xy(:,1)).*180./pi; % RSHO C7 ANG

walk_LSHO_xz_ang = atan2( walk_LSHO_xz(:,2), walk_LSHO_xz(:,1)).*180./pi; % LSHO C7 ANG
walk_RSHO_xz_ang = atan2( walk_RSHO_xz(:,2), walk_RSHO_xz(:,1)).*180./pi; % RSHO C7 ANG

walk_SHO_xy_ang = walk_LSHO_xy_ang-walk_RSHO_xy_ang; % SHO XY ANG

walk_SHO_xz_ang = walk_LSHO_xz_ang-walk_RSHO_xz_ang; % SHO XZ ANG

walk_head_ang = atan2(walk_head_v2(:,2),walk_head_v2(:,1))-atan2(walk_head_v1(:,2),walk_head_v1(:,1)); % HEAD ANG

walk_back_ang = atan2(walk_back_v2(:,2),walk_back_v2(:,1))-atan2(walk_back_v1(:,2),walk_back_v1(:,1)); %  BACK ANG

% Distance from Hip
walk_rwrst_dist = (walk_RHIPy - walk_RWRSTy);
walk_lwrst_dist = (walk_LHIPy - walk_LWRSTy);

disp('    Angles Calculated');

%% CREATE SLICE INDEX LIST

last_indx = length(Time_Step);
if min(lto_indx)<min(rto_indx)
    slice_foot = -1;
    indxes = lto_indx;
else
    slice_foot = 1;
    indxes = rto_indx;
end

%% WALKING ANGLE SLICES

struct_walk_back_angle = struct;
angles = struct;

for i=1:length(indxes)-1
    if i==1
        transition.SHO_XY_ang = walk_SHO_xy_ang(indxes(i):indxes(i+1));
        transition.SHO_XY_ang_mean = mean(transition.SHO_XY_ang);
        transition.SHO_XY_ang_std = std(transition.SHO_XY_ang);

        transition.SHO_XZ_ang = walk_SHO_xz_ang(indxes(i):indxes(i+1));
        transition.SHO_XZ_ang_mean = mean(transition.SHO_XZ_ang);
        transition.SHO_XZ_ang_std = std(transition.SHO_XZ_ang);

        transition.head_ang = walk_head_ang(indxes(i):indxes(i+1));
        transition.head_ang_mean = mean(transition.head_ang);
        transition.head_ang_std = std(transition.head_ang);

        transition.back_ang = walk_back_ang(indxes(i):indxes(i+1));
        transition.back_ang_mean = mean(transition.back_ang);
        transition.back_ang_std = std(transition.back_ang);

        transition.rwrst_dist = walk_rwrst_dist(indxes(i):indxes(i+1));
        transition.rwrst_dist_mean = mean(transition.rwrst_dist);
        transition.rwrst_dist_std = std(transition.rwrst_dist);

        transition.lwrst_dist = walk_lwrst_dist(indxes(i):indxes(i+1));
        transition.lwrst_dist_mean = mean(transition.lwrst_dist);
        transition.lwrst_dist_std = std(transition.lwrst_dist);
    else
        angles(i-1).SHO_XY_ang = walk_SHO_xy_ang(indxes(i):indxes(i+1));
        angles(i-1).SHO_XY_ang_mean = mean(angles(i-1).SHO_XY_ang);
        angles(i-1).SHO_XY_ang_std = std(angles(i-1).SHO_XY_ang);

        angles(i-1).SHO_XZ_ang = walk_SHO_xz_ang(indxes(i):indxes(i+1));
        angles(i-1).SHO_XZ_ang_mean = mean(angles(i-1).SHO_XZ_ang);
        angles(i-1).SHO_XZ_ang_std = std(angles(i-1).SHO_XZ_ang);

        angles(i-1).head_ang = walk_head_ang(indxes(i):indxes(i+1));
        angles(i-1).head_ang_mean = mean(angles(i-1).head_ang);
        angles(i-1).head_ang_std = std(angles(i-1).head_ang);

        angles(i-1).back_ang = walk_back_ang(indxes(i):indxes(i+1));
        angles(i-1).back_ang_mean = mean(angles(i-1).back_ang);
        angles(i-1).back_ang_std = std(angles(i-1).back_ang);

        angles(i-1).rwrst_dist = walk_rwrst_dist(indxes(i):indxes(i+1));
        angles(i-1).rwrst_dist_mean = mean(angles(i-1).rwrst_dist);
        angles(i-1).rwrst_dist_std = std(angles(i-1).rwrst_dist);

        angles(i-1).lwrst_dist = walk_lwrst_dist(indxes(i):indxes(i+1));
        angles(i-1).lwrst_dist_mean = mean(angles(i-1).lwrst_dist);
        angles(i-1).lwrst_dist_std = std(angles(i-1).lwrst_dist);
    
    end
end

struct_walk_back_angle.leg = slice_foot;
struct_walk_back_angle.steps = length(angles)+1; 

struct_walk_back_angle.transition_SHO_XY_ang_mean = mean([transition.SHO_XY_ang_mean]);
struct_walk_back_angle.transition_SHO_XY_ang_mean_std = std([transition.SHO_XY_ang_mean]);
struct_walk_back_angle.transition_SHO_XY_ang_std = mean([transition.SHO_XY_ang_std]);

struct_walk_back_angle.transition_SHO_XZ_ang_mean = mean([transition.SHO_XZ_ang_mean]);
struct_walk_back_angle.transition_SHO_XZ_ang_mean_std = std([transition.SHO_XZ_ang_mean]);
struct_walk_back_angle.transition_SHO_XZ_ang_std = mean([transition.SHO_XZ_ang_std]);

struct_walk_back_angle.transition_head_ang_mean = mean([transition.head_ang_mean]);
struct_walk_back_angle.transition_head_ang_mean_std = std([transition.head_ang_mean]);
struct_walk_back_angle.transition_head_ang_std = mean([transition.head_ang_std]);

struct_walk_back_angle.transition_back_ang_mean = mean([transition.back_ang_mean]);
struct_walk_back_angle.transition_back_ang_mean_std = std([transition.back_ang_mean]);
struct_walk_back_angle.transition_back_ang_std = mean([transition.back_ang_std]);

struct_walk_back_angle.transition_rwrst_dist_mean = mean([transition.rwrst_dist_mean]);
struct_walk_back_angle.transition_rwrst_dist_mean_std = std([transition.rwrst_dist_mean]);
struct_walk_back_angle.transition_rwrst_dist_std = mean([transition.rwrst_dist_std]);

struct_walk_back_angle.transition_lwrst_dist_mean = mean([transition.lwrst_dist_mean]);
struct_walk_back_angle.transition_lwrst_dist_mean_std = std([transition.lwrst_dist_mean]);
struct_walk_back_angle.transition_lwrst_dist_std = mean([transition.lwrst_dist_std]);

assignin('base','angles',angles);

if isempty(fieldnames(angles))
    struct_walk_back_angle.SHO_XY_ang_mean = NaN;
    struct_walk_back_angle.SHO_XY_ang_mean_std = NaN;
    struct_walk_back_angle.SHO_XY_ang_std = NaN;

    struct_walk_back_angle.SHO_XZ_ang_mean = NaN;
    struct_walk_back_angle.SHO_XZ_ang_mean_std = NaN;
    struct_walk_back_angle.SHO_XZ_ang_std = NaN;

    struct_walk_back_angle.head_ang_mean = NaN;
    struct_walk_back_angle.head_ang_mean_std = NaN;
    struct_walk_back_angle.head_ang_std = NaN;

    struct_walk_back_angle.back_ang_mean = NaN;
    struct_walk_back_angle.back_ang_mean_std = NaN;
    struct_walk_back_angle.back_ang_std = NaN;

    struct_walk_back_angle.rwrst_dist_mean = NaN;
    struct_walk_back_angle.rwrst_dist_mean_std = NaN;
    struct_walk_back_angle.rwrst_dist_std = NaN;

    struct_walk_back_angle.lwrst_dist_mean = NaN;
    struct_walk_back_angle.lwrst_dist_mean_std = NaN;
    struct_walk_back_angle.lwrst_dist_std = NaN;
else
    struct_walk_back_angle.SHO_XY_ang_mean = mean([angles.SHO_XY_ang_mean]);
    struct_walk_back_angle.SHO_XY_ang_mean_std = std([angles.SHO_XY_ang_mean]);
    struct_walk_back_angle.SHO_XY_ang_std = mean([angles.SHO_XY_ang_std]);

    struct_walk_back_angle.SHO_XZ_ang_mean = mean([angles.SHO_XZ_ang_mean]);
    struct_walk_back_angle.SHO_XZ_ang_mean_std = std([angles.SHO_XZ_ang_mean]);
    struct_walk_back_angle.SHO_XZ_ang_std = mean([angles.SHO_XZ_ang_std]);

    struct_walk_back_angle.head_ang_mean = mean([angles.head_ang_mean]);
    struct_walk_back_angle.head_ang_mean_std = std([angles.head_ang_mean]);
    struct_walk_back_angle.head_ang_std = mean([angles.head_ang_std]);

    struct_walk_back_angle.back_ang_mean = mean([angles.back_ang_mean]);
    struct_walk_back_angle.back_ang_mean_std = std([angles.back_ang_mean]);
    struct_walk_back_angle.back_ang_std = mean([angles.back_ang_std]);

    struct_walk_back_angle.rwrst_dist_mean = mean([angles.rwrst_dist_mean]);
    struct_walk_back_angle.rwrst_dist_mean_std = std([angles.rwrst_dist_mean]);
    struct_walk_back_angle.rwrst_dist_std = mean([angles.rwrst_dist_std]);

    struct_walk_back_angle.lwrst_dist_mean = mean([angles.lwrst_dist_mean]);
    struct_walk_back_angle.lwrst_dist_mean_std = std([angles.lwrst_dist_mean]);
    struct_walk_back_angle.lwrst_dist_std = mean([angles.lwrst_dist_std]);
end

disp('    Angles Sliced');


%% PLOTS

figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk'; % cVec = colour vector;
for k = 1:length(angles)+1
    
    step_num = strcat("Step No: ",num2str(k));
    
    subplot(2,1,1);
    if k==1
        plot(1:length(transition.SHO_XY_ang),transition.SHO_XY_ang,[cVec(k) '-'],'DisplayName',step_num);
    elseif isempty(fieldnames(angles))
        continue
    else
        plot(1:length(angles(k-1).SHO_XY_ang),angles(k-1).SHO_XY_ang,[cVec(k) '-'],'DisplayName',step_num);
    end
    hold on;
    title("SHO XY Angle");
    legend('-DynamicLegend');
    legend('show');
    
    subplot(2,1,2);
    if k==1
        plot(1:length(transition.SHO_XZ_ang),transition.SHO_XZ_ang,[cVec(k) '-'],'DisplayName',step_num);
    elseif isempty(fieldnames(angles))
        continue
    else
        plot(1:length(angles(k-1).SHO_XZ_ang),angles(k-1).SHO_XZ_ang,[cVec(k) '-'],'DisplayName',step_num);
    end
    hold on;
    title("SHO XZ Angle");
    legend('-DynamicLegend');
    legend('show');
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_walking_shoulder_ang.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk'; % cVec = colour vector;
for k = 1:length(angles)+1
    
    step_num = strcat("Step No: ",num2str(k));
    
    subplot(2,1,1);
    if k==1
        plot(1:length(transition.head_ang),transition.head_ang,[cVec(k) '-'],'DisplayName',step_num);
    elseif isempty(fieldnames(angles))
        continue
    else
        plot(1:length(angles(k-1).head_ang),angles(k-1).head_ang,[cVec(k) '-'],'DisplayName',step_num);
    end
    hold on;
    title("Head Angle");
    legend('-DynamicLegend');
    legend('show');
    
    subplot(2,1,2);
    if k==1
        plot(1:length(transition.back_ang),transition.back_ang,[cVec(k) '-'],'DisplayName',step_num);
    elseif isempty(fieldnames(angles))
        continue
    else
        plot(1:length(angles(k-1).back_ang),angles(k-1).back_ang,[cVec(k) '-'],'DisplayName',step_num);
    end
    hold on;
    title("Back Angle");
    legend('-DynamicLegend');
    legend('show');
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_walking_head_and_back_ang.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk'; % cVec = colour vector;
for k = 1:length(angles)+1
    
    step_num = strcat("Step No: ",num2str(k));
    
    subplot(2,1,1);
    if k==1
        plot(1:length(transition.rwrst_dist),transition.rwrst_dist,[cVec(k) '-'],'DisplayName',step_num);
    elseif isempty(fieldnames(angles))
        continue
    else
        plot(1:length(angles(k-1).rwrst_dist),angles(k-1).rwrst_dist,[cVec(k) '-'],'DisplayName',step_num);
    end
    hold on;
    title("Right Arm Distance");
    legend('-DynamicLegend');
    legend('show');
    
    subplot(2,1,2);
    if k==1
        plot(1:length(transition.lwrst_dist),transition.lwrst_dist,[cVec(k) '-'],'DisplayName',step_num);
    elseif isempty(fieldnames(angles))
        continue
    else
        plot(1:length(angles(k-1).lwrst_dist),angles(k-1).lwrst_dist,[cVec(k) '-'],'DisplayName',step_num);
    end
    hold on;
    title("Left Arm Distance");
    legend('-DynamicLegend');
    legend('show');
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_walking_dist_from_body.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp("    Angles Plotted");

%% SAVE .MAT FILE

if init_struct.mat_save
    
    save(strcat(p_num,'_walk_back_angle.mat'),'struct_walk_back_angle');
    disp('    .MAT Saved');

end

%% CHANGE BACK TO DIRECTORY

cd('..')

end