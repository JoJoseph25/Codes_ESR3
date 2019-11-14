function [gait_init_start,gait_init_end] = gait_init_calc(init_struct,p_num,ind_end,FP_data,Mrk_Data)
% function [gait_init_start,gait_init_end] = gait_init_detect(init_struct,p_num,Mrk_Data,FP_data)
% 
%  This function calculates the gait initialization (first step) from 
%  forceplate and marker data and outputs the start and end index
%  to slice the data and saves plot (if wanted).
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Number (string) 
%           Mrk_Data - Marker Data (table)
%           FP_data - Force plate data (table)
%     
%   OUTPUT: gait_ind_start - Start Index to slice data (index)
%           gait_ind_end - End Index to slice data (index)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% GAIT INITIALIZATION DETECTION

indx_frst_stp = ind_end + 600; % ind_end = 5 sec before first step

gait_init_start = indx_frst_stp - 120; % 1 sec before step

gait_init_end = indx_frst_stp + 600; % 5 sec after step

time_frst_stp = Mrk_Data.Time(indx_frst_stp);

%%

% Base of spine / Iliopelviic (IP)
mrkr_IPx = Mrk_Data.IPX(:); % IPX

% Right Ankle (RANK)
mrkr_RANKz = Mrk_Data.RANKZ(:); % RANKZ
mrkr_RANKx = Mrk_Data.RANKX(:); % RANKX

% Left Ankle (LANK)
mrkr_LANKz = Mrk_Data.LANKZ(:); % LANKZ
mrkr_LANKx = Mrk_Data.LANKX(:); % LANKX

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
mrkr_IPx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_IPx)/1000; % m

mrkr_RANKz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKz)/1000; % m
mrkr_RANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKx)/1000; % m

mrkr_LANKz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKz)/1000; % m
mrkr_LANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKx)/1000; % m

mrkr_RHEELz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELz)/1000; % m
mrkr_RHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELx)/1000; % m

mrkr_LHEELz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELz)/1000; % m
mrkr_LHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELx)/1000; % m

mrkr_RTOEz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEz)/1000; % m
mrkr_RTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEx)/1000; % m

mrkr_LTOEz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEz)/1000; % m
mrkr_LTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEx)/1000; % m

%% CLEAN MARKER DATA

% Base of spine / Iliopelviic (IP)
IPX = mrkr_IPx(gait_init_start:gait_init_end); % IPX

% Right Ankle (RANK)
RANKZ = mrkr_RANKz(gait_init_start:gait_init_end); % RANKZ
RANKX = mrkr_RANKx(gait_init_start:gait_init_end); % RANKX

% Left Ankle (LANK)
LANKZ = mrkr_LANKz(gait_init_start:gait_init_end); % LANKZ
LANKX = mrkr_LANKx(gait_init_start:gait_init_end); % LANKX

% Right Heel (RHEEL)
RHEELZ = mrkr_RHEELz(gait_init_start:gait_init_end); % RHEELZ
RHEELX = mrkr_RHEELx(gait_init_start:gait_init_end); % RHEELX

% Left Heel (LHEEL)
LHEELZ = mrkr_LHEELz(gait_init_start:gait_init_end); % LHEELZ
LHEELX = mrkr_LHEELx(gait_init_start:gait_init_end); % LHEELX

% Right Toe (RTOE)
RTOEZ = mrkr_RTOEz(gait_init_start:gait_init_end); % RTOEZ
RTOEX = mrkr_RTOEx(gait_init_start:gait_init_end); % RTOEX

% Left Toe (LTOE)
LTOEZ = mrkr_LTOEz(gait_init_start:gait_init_end); % LTOEZ
LTOEX = mrkr_LTOEx(gait_init_start:gait_init_end); % LTOEX

clear mrkr*

%% DETECT STEPS

RHeel_stps_dist = (IPX-RHEELX);
LHeel_stps_dist = (IPX-LHEELX);

RToe_stps_dist = (IPX-RTOEX);
LToe_stps_dist = (IPX-LTOEX);

[peak_rstp, peak_rstp_indx] = findpeaks(RHeel_stps_dist,'MinPeakDistance',100);
[peak_lstp, peak_lstp_indx] = findpeaks(LHeel_stps_dist,'MinPeakDistance',100);

[valley_rstp, valley_rstp_indx] = findpeaks(RToe_stps_dist,'MinPeakDistance',100);
[valley_lstp, valley_lstp_indx] = findpeaks(LToe_stps_dist,'MinPeakDistance',100);


% Time vector
Time_Step = Mrk_Data.Time(gait_init_start:gait_init_end);
find(Time_Step==time_frst_stp,1,'first')

%% DISTANCE

figure;
plot(Mrk_Data.Time(gait_init_start:gait_init_end),RHeel_stps_dist);
hold on;
line([time_frst_stp time_frst_stp],[-500 500],'Color',[1 0 0]);
plot(Time_Step(peak_rstp_indx),RHeel_stps_dist(peak_rstp_indx),'*');
hold off;
title("IP to RHEEL ");
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'IP_2_RHeel.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end


figure;
plot(Mrk_Data.Time(gait_init_start:gait_init_end),LHeel_stps_dist);
hold on;
line([time_frst_stp time_frst_stp],[-500 500],'Color',[1 0 0]);
plot(Time_Step(peak_lstp_indx),LHeel_stps_dist(peak_lstp_indx),'*');
hold off;
title("IP to LHEEL ");
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'IP_2_LHeel.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%%
peak_rstp_indx=peak_rstp_indx(1:4);
peak_lstp_indx=peak_lstp_indx(1:4);
%% HEEL

figure;
subplot(3,1,1);
plot(FP_data.TIME(gait_init_start:gait_init_end),FP_data.Force_Z(gait_init_start:gait_init_end));
line([time_frst_stp time_frst_stp],[-1000 1000],'Color',[1 0 0])
title("Force Z across Recording")
ylabel('Force Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,2);
plot(Mrk_Data.Time(gait_init_start:gait_init_end),RHEELX);
hold on;
line([time_frst_stp time_frst_stp],[-100 2000],'Color',[1 0 0]);
plot(Time_Step(peak_rstp_indx),RHEELX(peak_rstp_indx),'*');
hold off;
title("RHEEL across Recording");
ylabel('RHEEL X-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,3);
plot(Mrk_Data.Time(gait_init_start:gait_init_end),LHEELX);
hold on;
line([time_frst_stp time_frst_stp],[-100 2000],'Color',[1 0 0]);
plot(Time_Step(peak_lstp_indx),LHEELX(peak_lstp_indx),'*');
hold off;
title("LHEEL across Recording");
ylabel('LHEEL X-Axis'),xlabel(strcat('Time (secs)'));
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'Heel_gait_init.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end


%% HEEL X-Z

figure;
subplot(2,1,1);
plot(RHEELX,RHEELZ);
hold on;
plot(RHEELX(peak_rstp_indx),RHEELZ(peak_rstp_indx),'*');
hold off;
title("RHEEL X vs Z");
ylabel('RHEEL X-Axis'),xlabel('RHEEL X-Axis');

subplot(2,1,2);
plot(LHEELX,LHEELZ);
hold on;
plot(LHEELX(peak_lstp_indx),LHEELZ(peak_lstp_indx),'*');
hold off;
title("LHEEL X vs Z");
ylabel('LHEEL Z-Axis'),xlabel('LHEEL X-Axis');


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'Heel_XZ_gait_init.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% ANKLE

% RANK
figure;
subplot(3,1,1);
plot(FP_data.TIME(gait_init_start:gait_init_end),FP_data.Force_Z(gait_init_start:gait_init_end));
line([time_frst_stp time_frst_stp],[-1000 1000],'Color',[1 0 0])
title("Force Z across Recording")
ylabel('Force Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,2);
plot(Mrk_Data.Time(gait_init_start:gait_init_end),RANKZ);
hold on;
line([time_frst_stp time_frst_stp],[0 500],'Color',[1 0 0]);
plot(Time_Step(peak_rstp_indx),RANKZ(peak_rstp_indx),'*');
hold off;
title("RANK across Recording");
ylabel('RANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,3);
plot(Mrk_Data.Time(gait_init_start:gait_init_end),RANKX);
hold on;
line([time_frst_stp time_frst_stp],[-2000 2000],'Color',[1 0 0]);
plot(Time_Step(peak_rstp_indx),RANKX(peak_rstp_indx),'*');
hold off;
title("RANK across Recording");
ylabel('RANK X-Axis'),xlabel(strcat('Time (secs)'));

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'RANK_gait_init.jpg'));
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
plot(Mrk_Data.Time(gait_init_start:gait_init_end),LANKZ);
hold on;
line([time_frst_stp time_frst_stp],[0 500],'Color',[1 0 0]);
plot(Time_Step(peak_lstp_indx),LANKZ(peak_lstp_indx),'*');
hold off;
title("LANK across Recording");
ylabel('LANK Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,3);
plot(Mrk_Data.Time(gait_init_start:gait_init_end),LANKX);
hold on;
line([time_frst_stp time_frst_stp],[-2000 2000],'Color',[1 0 0]);
plot(Time_Step(peak_lstp_indx),LANKX(peak_lstp_indx),'*');
hold off;
title("LANK across Recording");
ylabel('LANK X-Axis'),xlabel(strcat('Time (secs)'));


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'LANK_gait_init.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% ANKLE X-Z

figure;
subplot(2,1,1);
plot(RANKX,RANKZ);
hold on;
plot(RANKX(peak_rstp_indx),RANKZ(peak_rstp_indx),'*');
hold off;
title("RANK X vs Z");
ylabel('RANK X-Axis'),xlabel('RANK X-Axis');

subplot(2,1,2);
plot(LANKX,LANKZ);
hold on;
plot(LANKX(peak_lstp_indx),LANKZ(peak_lstp_indx),'*');
hold off;
title("LANK X vs Z");
ylabel('LANK Z-Axis'),xlabel('LANK X-Axis');


if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'ANK_XZ_gait_init.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

end