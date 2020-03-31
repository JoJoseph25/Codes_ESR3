function [ind_start,ind_end] = first_step_detect(init_struct,p_num,Mrk_Data,FP_data)
% function [ind_start,ind_end] = first_step_detect(p_num,Mrk_Data,FP_data)
% 
%  This function calculates the first step from 
%  forceplate and marker data and outputs the start and end index
%  to slice the data and saves plot (if wanted).
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Number (string) 
%           Mrk_Data - Marker Data (table)
%           FP_data - Force plate data (table)
%     
%   OUTPUT: ind_start - Start Index to slice data (index)
%           ind_end - End Index to slice data (index)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% EXCEPTIONS

exceptions ={'1020_dy5';};

%% TIME
f10_sec = 10*20;

f30_sec=30*120; % first 30-sec data point (index) 

mrk_time= Mrk_Data.Time(:); % Marker Data Time

fp_time = FP_data.TIME(:); % Forceplate Data Time

%% DETECTION

%Data to be used for detection
Force_Z = FP_data.Force_Z(f30_sec:end); % Force Z
RHEEL_X = Mrk_Data.RHEELX(f30_sec:end); % RHEEL X
LHEEL_X = Mrk_Data.LHEELX(f30_sec:end); % LHEEL X

if ismember(p_num,exceptions)
    detect_criteria = 465;
else
    detect_criteria = -5;
end

% Value of Force Z around zero (first value Force_Z>+/-5 + RHEEL\LHEEL>50mm)
indx_frst_stp = find((Force_Z > detect_criteria & (RHEEL_X > 50 | LHEEL_X > 50)),1,'first'); % Index First Step of Force Plate

indx_frst_stp = indx_frst_stp + f30_sec;     

%% TIME POINT

% Time point when the subject moves from Force Plate
time_frst_stp = fp_time(indx_frst_stp); % Time First Step Initiation

% Time Range
r1 = 30; % Start Time for data slice
r2 = fp_time(indx_frst_stp)-5; % Time range end point [sec]
t_rng = r2 - r1; % Everything will be calculated for time range

% Additional Indicators 
min_time = min(fp_time); % Start Time
max_time = max(fp_time); % End Time
t_len = max_time - min_time; % Time Length 

ssfs = time_frst_stp-min_time; % Seconds from Start to First Step 
sfse = max_time-time_frst_stp; % Seconds from First Step till E

%% INDEXES

% Index
indx_bfr_frst_stp = 3600; % 30 seconds in terms of index i.e. 30*120 = 3600
indx_frst_stp_strt = indx_frst_stp-indx_bfr_frst_stp;  % 30 seconds before first step
indx_frst_stp_end = indx_frst_stp-600; % 5 sec before first step

% Define the Start
ind_start = f30_sec; % data after first 30 seconds

% Define the End
ind_end=indx_frst_stp_end; %5 sec before first step initiation

%% PLOT FIRST STEP

figure;
subplot(3,1,1);
plot(FP_data.TIME,FP_data.Force_Z);
line([time_frst_stp time_frst_stp],[-1000 1000],'Color',[1 0 0])
title("Force Z across Recording")
ylabel('Force Z-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,2);
plot(Mrk_Data.Time,Mrk_Data.RHEELX);
line([time_frst_stp time_frst_stp],[-5000 5000],'Color',[1 0 0])
title("RHEEL across Recording")
ylabel('RHEEL X-Axis'),xlabel(strcat('Time (secs)'));

subplot(3,1,3);
plot(Mrk_Data.Time,Mrk_Data.LHEELX);
line([time_frst_stp time_frst_stp],[-5000 5000],'Color',[1 0 0])
title("LHEEL across Recording")
ylabel('LHEEL X-Axis'),xlabel(strcat('Time (secs)'));

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_first_step_detect.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% DISPLAY

disp('First Step Detection'); % print out to command line

end