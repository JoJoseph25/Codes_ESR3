function [struct_stand_angle] = stand_angle(init_struct,p_num,Mrk_Data,ind_start,ind_end)
% function [struct_stand_angle] = stand_angle(init_struct,p_num,FP_data,ind_start,ind_end)
% 
%  This function calculates the body angles while saving
%  struct_stand_angle and plots all the graph and .mat file
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Condition Number (string)
%           Mrk_Data - Marker data (table)
%           ind_start - Start Index to slice data (index)
%           ind_end - End Index to slice data (index)
%     
%   OUTPUT: struct_stand_COP - COP parameters (struct)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% EXCEPTIONS

% exceptions={'1010_dy6','1010_dy5'}; % List of exceptions that need different 
%                                     min distance than normal data for
%                                     clustering      

%% CHANGE FOLDER

if ~exist('Stand Angle', 'dir') % Check if folder exist
    mkdir('Stand Angle'); % make new folder
end

cd('Stand Angle'); % Change directory to new folder

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

% Marker Data Standing
stand_range=num2str((ind_end-ind_start)/init_struct.Fs);% Time range 

stand_time = Mrk_Data.Time(ind_start:ind_end,:); % Time vector

% Top of Head (TPHD)
stand_TPHDx = mrkr_TPHDx_filter(ind_start:ind_end);
stand_TPHDy = mrkr_TPHDy_filter(ind_start:ind_end);
stand_TPHDz = mrkr_TPHDz_filter(ind_start:ind_end);

% Base of Neck / Cervical 7 (C7)
stand_C7x = mrkr_C7x_filter(ind_start:ind_end);
stand_C7y = mrkr_C7y_filter(ind_start:ind_end);
stand_C7z = mrkr_C7z_filter(ind_start:ind_end);

% Right Shoulder (RSHO)
stand_RSHOx = mrkr_RSHOx_filter(ind_start:ind_end);
stand_RSHOy = mrkr_RSHOy_filter(ind_start:ind_end);
stand_RSHOz = mrkr_RSHOz_filter(ind_start:ind_end);

% Left Shoulder (SHO)
stand_LSHOx = mrkr_LSHOx_filter(ind_start:ind_end);
stand_LSHOy = mrkr_LSHOy_filter(ind_start:ind_end);
stand_LSHOz = mrkr_LSHOz_filter(ind_start:ind_end);

% Base of spine / Iliopelviic (IP)
stand_IPx = mrkr_IPx_filter(ind_start:ind_end);
stand_IPy = mrkr_IPy_filter(ind_start:ind_end);
stand_IPz = mrkr_IPz_filter(ind_start:ind_end);

% Right Wrist (RWRST)
stand_RWRSTz = mrkr_RWRSTz_filter(ind_start:ind_end);
stand_RWRSTy = mrkr_RWRSTy_filter(ind_start:ind_end);

% Left Wrist (LWRST)
stand_LWRSTz = mrkr_LWRSTz_filter(ind_start:ind_end);
stand_LWRSTy = mrkr_LWRSTy_filter(ind_start:ind_end);

% Left HIP (LHIP)
stand_LHIPy = mrkr_LHIPy_filter(ind_start:ind_end);

% Right HIP (LHIP)
stand_RHIPy = mrkr_RHIPy_filter(ind_start:ind_end);

clear mrkr_* % clear unnecessary variables 

disp('Angles: ');

%% Angle Calculation
stand_RSHO_xy = [stand_RSHOx-stand_C7x stand_RSHOy-stand_C7y];
stand_LSHO_xy = [stand_LSHOx-stand_C7x stand_LSHOy-stand_C7y];

stand_RSHO_xz = [stand_RSHOx-stand_C7x stand_RSHOz-stand_C7z];
stand_LSHO_xz = [stand_LSHOx-stand_C7x stand_LSHOz-stand_C7z];

stand_head_v1 = [stand_TPHDx-stand_C7x stand_TPHDz-stand_C7z];
stand_head_v2 = [(stand_TPHDx-stand_C7x)*0 stand_TPHDz-stand_C7z];

stand_back_v1 = [stand_IPx-stand_C7x stand_IPz-stand_C7z];
stand_back_v2 = [(stand_IPx-stand_C7x)*0 stand_IPz-stand_C7z];

% Upper Body Angles
stand_LSHO_xy_ang = atan2( stand_LSHO_xy(:,2), stand_LSHO_xy(:,1)).*180./pi; % LSHO C7 ANG
stand_RSHO_xy_ang = atan2( stand_RSHO_xy(:,2), stand_RSHO_xy(:,1)).*180./pi; % RSHO C7 ANG
 
% neg_LSHO_ang_indx = find(stand_LSHO_xy_ang<0);
% 
% stand_LSHO_ang(neg_LSHO_ang_indx)

stand_LSHO_xz_ang = atan2( stand_LSHO_xz(:,2), stand_LSHO_xz(:,1)).*180./pi; % LSHO C7 ANG
stand_RSHO_xz_ang = atan2( stand_RSHO_xz(:,2), stand_RSHO_xz(:,1)).*180./pi; % RSHO C7 ANG

% Replace -ve angles with 360 + (-ve ang)
stand_LSHO_xy_ang(stand_LSHO_xy_ang<0) = 360.+ stand_LSHO_xy_ang(stand_LSHO_xy_ang<0); 
stand_LSHO_xz_ang(stand_LSHO_xz_ang<0) = 360.+ stand_LSHO_xz_ang(stand_LSHO_xz_ang<0);

stand_RSHO_xy_ang(stand_RSHO_xy_ang<0) = 360.+ stand_RSHO_xy_ang(stand_RSHO_xy_ang<0); 
stand_RSHO_xz_ang(stand_RSHO_xz_ang<0) = 360.+ stand_RSHO_xz_ang(stand_RSHO_xz_ang<0);

stand_SHO_xy_ang = stand_LSHO_xy_ang-stand_RSHO_xy_ang; % SHO XY ANG
stand_SHO_xz_ang = stand_LSHO_xz_ang-stand_RSHO_xz_ang; % SHO XZ ANG

stand_head_ang = atan2(stand_head_v2(:,2),stand_head_v2(:,1))-atan2(stand_head_v1(:,2),stand_head_v1(:,1)); % HEAD ANG

stand_back_ang = atan2(stand_back_v2(:,2),stand_back_v2(:,1))-atan2(stand_back_v1(:,2),stand_back_v1(:,1)); %  BACK ANG

% Distance from Hip
rwrst_dist = (stand_RHIPy - stand_RWRSTy);
lwrst_dist = (stand_LHIPy - stand_LWRSTy);

disp('    Calculations Done');

%% CLUSTERING
% This is to cluster when you raise hand and keep portion of stable standing

rwrst_clust_mat = horzcat(stand_RWRSTy,stand_RWRSTz); % right wrist matrix
lwrst_clust_mat = horzcat(stand_LWRSTy,stand_LWRSTz); % left wrist matrix

rwrst_clust_mat = transpose(rwrst_clust_mat); % right wrist matrix transpose
lwrst_clust_mat = transpose(lwrst_clust_mat); % left wrist matrix transpose

assignin('base','rwrst_clust_mat',rwrst_clust_mat);
assignin('base','lwrst_clust_mat',lwrst_clust_mat);

min_dist= .18; % min distance for cluster

% Function of Mean Shift Cluster to form initial clusters on entire standing Right Wrist data
[rwrst_clustCent,rwrst_point2cluster,rwrst_clustMembsCell] = MeanShiftCluster(rwrst_clust_mat,min_dist);
% Function of Mean Shift Cluster to form initial clusters on entire standing Left Wrist data
[lwrst_clustCent,lwrst_point2cluster,lwrst_clustMembsCell] = MeanShiftCluster(lwrst_clust_mat,min_dist);
rwrst_numClust = length(rwrst_clustMembsCell); % number of clusters left wrist
lwrst_numClust = length(lwrst_clustMembsCell); % number of clusters right wrist

%Right Wrist Cleaning
for k=1:rwrst_numClust
    if length(rwrst_clustMembsCell{k})>240
        rwrst_indx_temp = rwrst_clustMembsCell{k}; % current right wrist cluster indexes
        rwrst_dist_temp = rwrst_clust_mat(:,rwrst_clustMembsCell{k}); % current right wrist cluster data points 
        rwrst_dist_temp = transpose(rwrst_dist_temp); % transpose the current right wrist cluster data points

        % create distance vector to be used to clean data or remove transitions
        rwrst_dist_mat = pdist(rwrst_dist_temp,'mahalanobis'); % Malabonis distance vector 
        clear rwrst_dist_temp % remove temp data point matrix as it is not needed
        rwrst_dist_mat = squareform(rwrst_dist_mat); % symmetric rwrist distance matrix

        Eps = 0.65; % min distance 
        MinPts = 240; % min data points need to call as cluster (1 sec=120*1 points)

        % Function of Density based clustering to remove the transitions 
        rwrst_Clust = DBSCAN(rwrst_dist_mat,Eps,MinPts);
        rwrst_indx_clust{k,1} = rwrst_indx_temp(rwrst_Clust~=0); % indexes of clean clusters    
    end
end

%Left Wrist Cleaning
for k=1:lwrst_numClust
    if length(lwrst_clustMembsCell{k})>240
        lwrst_indx_temp = lwrst_clustMembsCell{k}; % current left wrist cluster indexes
        lwrst_dist_temp = lwrst_clust_mat(:,lwrst_clustMembsCell{k}); % current left wrist cluster data points 
        lwrst_dist_temp = transpose(lwrst_dist_temp); % transpose the current left wrist cluster data points

        % create distance vector to be used to clean data or remove transitions
        lwrst_dist_mat = pdist(lwrst_dist_temp,'mahalanobis'); % Malabonis distance vector
        clear lwrst_dist_temp % remove temp data point matrix as it is not needed
        lwrst_dist_mat = squareform(lwrst_dist_mat); % symmetric rwrist distance matrix

        Eps = 0.65; % min distance 
        MinPts = 240; % min data points need to call as cluster (1 sec=120*1 points)

        % Function of Density based clustering to remove the transitions 
        lwrst_Clust = DBSCAN(lwrst_dist_mat,Eps,MinPts);
        lwrst_indx_clust{k,1} = lwrst_indx_temp(lwrst_Clust~=0); % indexes of clean clusters
    end
end

% Choose smallest cluster and keep the other for finding points common to both
if length(lwrst_indx_clust)>length(rwrst_indx_clust)
    
    min_wrst_unclean = rwrst_clustMembsCell; % Right wrist unclean
    min_wrst_indx_clust = rwrst_indx_clust; % Right wrist minimum 
    comp_indx_clust = horzcat(lwrst_indx_clust{:}); % Left wrist clean indexes

elseif length(rwrst_indx_clust)>length(lwrst_indx_clust)
    
    min_wrst_unclean = lwrst_clustMembsCell; % Left wrist unclean
    min_wrst_indx_clust = lwrst_indx_clust; % Left wrist minimum
    comp_indx_clust = horzcat(rwrst_indx_clust{:}); % Right wrist clean indexes

else
    min_wrst_unclean = rwrst_clustMembsCell; % Right wrist unclean
    min_wrst_indx_clust = rwrst_indx_clust; % Right wrist minimum
    comp_indx_clust = horzcat(lwrst_indx_clust{:}); % Left wrist clean indexes
end

% Keep index common to both right and left wrist clusters
for i=1:length(min_wrst_indx_clust)
    clust_wrst_indx{i,1} = intersect(min_wrst_indx_clust{i,1},comp_indx_clust);
end

full_wrst_indx = horzcat(clust_wrst_indx{:}); % Clean Wrist combined indexes

disp('    Clustering');

%% PLOTS

% Unclean
figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk'; % cVec = colour vector;
for k = 1:length(min_wrst_indx_clust)
    
    wrst_num=strcat("Wrst Clust ",num2str(k));
        
    subplot(2,1,1);
    hold on;
    plot(stand_time(min_wrst_unclean{k}),stand_RWRSTz(min_wrst_unclean{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("RIGHT WRIST ")
    ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;

    subplot(2,1,2);
    hold on;
    plot(stand_time(min_wrst_unclean{k}),stand_LWRSTz(min_wrst_unclean{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("LEFT WRIST")
    ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_wrist_elevation.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Clean
figure;
for k = 1:length(clust_wrst_indx)
    
    wrst_num=strcat("Wrst Clust ",num2str(k));
    
    subplot(2,1,1);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),stand_RWRSTz(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("RIGHT WRIST Clean")
    ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;

    subplot(2,1,2);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),stand_LWRSTz(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("LEFT WRIST Clean")
    ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'clean_wrist_elevation.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Wrist Plots');

%% SHOULDER ANGLES

figure;
for k = 1:length(clust_wrst_indx)
    
    wrst_num=strcat("Wrst Clust ",num2str(k));
    
    subplot(2,1,1);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),stand_SHO_xy_ang(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("Shoulder XY Angle")
    ylabel('Shoulder XY'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;

    subplot(2,1,2);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),stand_SHO_xz_ang(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("Shoulder XZ Angle")
    ylabel('Shoulder XZ'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;
    
end    
    
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'shoulder_ang.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% NECK & BACK ANGLES

figure;
for k = 1:length(clust_wrst_indx)
    
    wrst_num=strcat("Wrst Clust ",num2str(k));
    subplot(2,1,1);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),stand_head_ang(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("Neck Angle")
    ylabel('Neck'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;

    subplot(2,1,2);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),stand_back_ang(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("Back Angle")
    ylabel('Back'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'back_and_neck_ang.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% DISTANCE FROM BODY

figure;
for k = 1:length(clust_wrst_indx)
    
    wrst_num=strcat("Wrst Clust ",num2str(k));
    subplot(2,1,1);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),rwrst_dist(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("Right Hand Distance From Body")
    ylabel('RHand dist from body'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;

    subplot(2,1,2);
    hold on;
    plot(stand_time(clust_wrst_indx{k}),lwrst_dist(clust_wrst_indx{k}),[cVec(k) '.'],'DisplayName',wrst_num);
    title("Left Hand Distance From Body")
    ylabel('LHand dist from body'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'dist_from_body.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% ANGLE PARAMETERS 

clust_struct_stand_angle = struct;
struct_stand_angle = struct;

if length(clust_wrst_indx)>1 % More than one cluster
    for k=1:length(clust_wrst_indx)
        if length(clust_wrst_indx{k})>1
            %% Cluster Angle Parameters
            
            % Cluster info
            clust_struct_stand_angle(k).indx = clust_wrst_indx{k}'; % index
            clust_struct_stand_angle(k).data_pnts = length(clust_wrst_indx{k}); % data points
            clust_struct_stand_angle(k).time = length(clust_wrst_indx{k})/init_struct.Fs; % time
            
            % Angles
            clust_struct_stand_angle(k).clust_SHO_xy_ang = stand_SHO_xy_ang(clust_struct_stand_angle(k).indx); % SHOULDER XY angle of specific cluster in degrees
            clust_struct_stand_angle(k).clust_SHO_xz_ang = stand_SHO_xz_ang(clust_struct_stand_angle(k).indx); % SHOULDER XZ angle of specific cluster in degrees
            clust_struct_stand_angle(k).clust_head_ang = stand_head_ang(clust_struct_stand_angle(k).indx); % HEAD angle of specific cluster in degrees
            clust_struct_stand_angle(k).clust_back_ang = stand_back_ang(clust_struct_stand_angle(k).indx); % BACK angle of specific cluster in degrees
            clust_struct_stand_angle(k).clust_rwrst_dist = rwrst_dist(clust_struct_stand_angle(k).indx); % RHAND dist of specific cluster in meters
            clust_struct_stand_angle(k).clust_lwrst_dist = lwrst_dist(clust_struct_stand_angle(k).indx); % RHAND dist of specific cluster in meters
            
            % Mean
            clust_struct_stand_angle(k).clust_mean_SHO_xy_ang = mean(clust_struct_stand_angle(k).clust_SHO_xy_ang);
            clust_struct_stand_angle(k).clust_mean_SHO_xz_ang = mean(clust_struct_stand_angle(k).clust_SHO_xz_ang);
            clust_struct_stand_angle(k).clust_mean_head_ang = mean(clust_struct_stand_angle(k).clust_head_ang);
            clust_struct_stand_angle(k).clust_mean_back_ang = mean(clust_struct_stand_angle(k).clust_back_ang);
            clust_struct_stand_angle(k).clust_mean_rwrst_dist = mean(clust_struct_stand_angle(k).clust_rwrst_dist);
            clust_struct_stand_angle(k).clust_mean_lwrst_dist = mean(clust_struct_stand_angle(k).clust_lwrst_dist);
            
            % Standard-deviation
            clust_struct_stand_angle(k).clust_std_SHO_xy_ang = std(clust_struct_stand_angle(k).clust_SHO_xy_ang);
            clust_struct_stand_angle(k).clust_std_SHO_xz_ang = std(clust_struct_stand_angle(k).clust_SHO_xz_ang);
            clust_struct_stand_angle(k).clust_std_head_ang = std(clust_struct_stand_angle(k).clust_head_ang);
            clust_struct_stand_angle(k).clust_std_back_ang = std(clust_struct_stand_angle(k).clust_back_ang);
            clust_struct_stand_angle(k).clust_std_rwrst_dist = std(clust_struct_stand_angle(k).clust_rwrst_dist);
            clust_struct_stand_angle(k).clust_std_lwrst_dist = std(clust_struct_stand_angle(k).clust_lwrst_dist);
                        
        end
    end 
% Only one cluster     
else
    % Cluster info
    clust_struct_stand_angle.indx = full_wrst_indx'; % index
    clust_struct_stand_angle.data_pnts = length(full_wrst_indx); % data points
    clust_struct_stand_angle.time = length(full_wrst_indx)/init_struct.Fs; % time
    
    % Angles
    clust_struct_stand_angle.clust_SHO_xy_ang = stand_SHO_xy_ang(clust_struct_stand_angle.indx); % SHOULDER XY angle of specific cluster in degrees
    clust_struct_stand_angle.clust_SHO_xz_ang = stand_SHO_xz_ang(clust_struct_stand_angle.indx); % SHOULDER XZ angle of specific cluster in degrees
    clust_struct_stand_angle.clust_head_ang = stand_head_ang(clust_struct_stand_angle.indx); % HEAD angle of specific cluster in degrees
    clust_struct_stand_angle.clust_back_ang = stand_back_ang(clust_struct_stand_angle.indx); % BACK angle of specific cluster in degrees
    clust_struct_stand_angle.clust_rwrst_dist = rwrst_dist(clust_struct_stand_angle.indx); % RHAND dist of specific cluster in meters
    clust_struct_stand_angle.clust_lwrst_dist = lwrst_dist(clust_struct_stand_angle.indx); % RHAND dist of specific cluster in meters

    % Mean
    clust_struct_stand_angle.clust_mean_SHO_xy_ang = mean(clust_struct_stand_angle.clust_SHO_xy_ang);
    clust_struct_stand_angle.clust_mean_SHO_xz_ang = mean(clust_struct_stand_angle.clust_SHO_xz_ang);
    clust_struct_stand_angle.clust_mean_head_ang = mean(clust_struct_stand_angle.clust_head_ang);
    clust_struct_stand_angle.clust_mean_back_ang = mean(clust_struct_stand_angle.clust_back_ang);
    clust_struct_stand_angle.clust_mean_rwrst_dist = mean(clust_struct_stand_angle.clust_rwrst_dist);
    clust_struct_stand_angle.clust_mean_lwrst_dist = mean(clust_struct_stand_angle.clust_lwrst_dist);

    % Standard-deviation
    clust_struct_stand_angle.clust_std_SHO_xy_ang = std(clust_struct_stand_angle.clust_SHO_xy_ang);
    clust_struct_stand_angle.clust_std_SHO_xz_ang = std(clust_struct_stand_angle.clust_SHO_xz_ang);
    clust_struct_stand_angle.clust_std_head_ang = std(clust_struct_stand_angle.clust_head_ang);
    clust_struct_stand_angle.clust_std_back_ang = std(clust_struct_stand_angle.clust_back_ang);
    clust_struct_stand_angle.clust_std_rwrst_dist = std(clust_struct_stand_angle.clust_rwrst_dist);
    clust_struct_stand_angle.clust_std_lwrst_dist = std(clust_struct_stand_angle.clust_lwrst_dist);
                        
end

%% Full COP Parameters

struct_stand_angle.full_indx = full_wrst_indx';
struct_stand_angle.time = length(full_wrst_indx)/init_struct.Fs; 

% Angles
struct_stand_angle.full_SHO_xy_ang = stand_SHO_xy_ang(full_wrst_indx); % SHOULDER XY angle in degrees
struct_stand_angle.full_SHO_xz_ang = stand_SHO_xz_ang(full_wrst_indx); % SHOULDER XZ angle in degrees
struct_stand_angle.full_head_ang = stand_head_ang(full_wrst_indx); % HEAD angle in degrees
struct_stand_angle.full_back_ang = stand_back_ang(full_wrst_indx); % BACK angle in degrees
struct_stand_angle.full_rwrst_dist = rwrst_dist(full_wrst_indx); % RHAND dist in meters
struct_stand_angle.full_lwrst_dist = lwrst_dist(full_wrst_indx); % RHAND dist in meters

 % Mean
struct_stand_angle.full_mean_SHO_xy_ang = mean(struct_stand_angle.full_SHO_xy_ang);
struct_stand_angle.full_mean_SHO_xz_ang = mean(struct_stand_angle.full_SHO_xz_ang);
struct_stand_angle.full_mean_head_ang = mean(struct_stand_angle.full_head_ang);
struct_stand_angle.full_mean_back_ang = mean(struct_stand_angle.full_back_ang);
struct_stand_angle.full_mean_rwrst_dist = mean(struct_stand_angle.full_rwrst_dist);
struct_stand_angle.full_mean_lwrst_dist = mean(struct_stand_angle.full_lwrst_dist);

% Standard-deviation
struct_stand_angle.full_std_SHO_xy_ang = std(struct_stand_angle.full_SHO_xy_ang);
struct_stand_angle.full_std_SHO_xz_ang = std(struct_stand_angle.full_SHO_xz_ang);
struct_stand_angle.full_std_head_ang = std(struct_stand_angle.full_head_ang);
struct_stand_angle.full_std_back_ang = std(struct_stand_angle.full_back_ang);
struct_stand_angle.full_std_rwrst_dist = std(struct_stand_angle.full_rwrst_dist);
struct_stand_angle.full_std_lwrst_dist = std(struct_stand_angle.full_lwrst_dist);

%% RELATIVE PARAMETERS

% parameters adjusted by cluster-size (n*P/n) where n= num of points and
%                                                   P= parameter of interest

% Relative Mean
struct_stand_angle.clust_rel_mean_SHO_xy_ang = (sum([clust_struct_stand_angle.clust_mean_SHO_xy_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_mean_SHO_xz_ang = (sum([clust_struct_stand_angle.clust_mean_SHO_xz_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_mean_head_ang = (sum([clust_struct_stand_angle.clust_mean_head_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_mean_back_ang = (sum([clust_struct_stand_angle.clust_mean_back_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_mean_rwrst_dist = (sum([clust_struct_stand_angle.clust_mean_rwrst_dist].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_mean_lwrst_dist = (sum([clust_struct_stand_angle.clust_mean_lwrst_dist].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);

% Relative Standard-deviation
struct_stand_angle.clust_rel_std_SHO_xy_ang = (sum([clust_struct_stand_angle.clust_std_SHO_xy_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_std_SHO_xz_ang = (sum([clust_struct_stand_angle.clust_std_SHO_xz_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_std_head_ang = (sum([clust_struct_stand_angle.clust_std_head_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_std_back_ang = (sum([clust_struct_stand_angle.clust_std_back_ang].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_std_rwrst_dist = (sum([clust_struct_stand_angle.clust_std_rwrst_dist].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);
struct_stand_angle.clust_rel_std_lwrst_dist = (sum([clust_struct_stand_angle.clust_std_lwrst_dist].*[clust_struct_stand_angle.data_pnts]))/sum([clust_struct_stand_angle.data_pnts]);

disp("    Parameters Calculated");

%% CONTINOUS PARAMETERS

% For clusters
clust_continous_angle=struct;

for i=1:length(clust_struct_stand_angle)
   clear cont*
   
    cont_SHOxy_ang = [clust_struct_stand_angle(i).clust_SHO_xy_ang,clust_struct_stand_angle(i).indx]; % create matrix of value and index
    cont_SHOxy_ang = sortrows(cont_SHOxy_ang,2); % sort according to values
    cont_SHO_xy_ang = rollstat(cont_SHOxy_ang,600,600); % calculate rolling window stats (120Hz * 5 sec = 600 points)
    clust_continous_angle(i).SHO_xy_ang = [cont_SHO_xy_ang.stats]; % save to struct

    cont_SHOxz_ang = [clust_struct_stand_angle(i).clust_SHO_xz_ang,clust_struct_stand_angle(i).indx];
    cont_SHOxz_ang = sortrows(cont_SHOxz_ang,2); 
    cont_SHO_xz_ang = rollstat(cont_SHOxz_ang,600,600);
    clust_continous_angle(i).SHO_xz_ang = [cont_SHO_xz_ang.stats];

    cont_Head_ang = [clust_struct_stand_angle(i).clust_head_ang,clust_struct_stand_angle(i).indx];
    cont_Head_ang = sortrows(cont_Head_ang,2); 
    cont_head_ang = rollstat(cont_Head_ang,600,600);
    clust_continous_angle(i).head_ang = [cont_head_ang.stats];

    cont_Back_ang = [clust_struct_stand_angle(i).clust_back_ang,clust_struct_stand_angle(i).indx];
    cont_Back_ang = sortrows(cont_Back_ang,2); 
    cont_back_ang = rollstat(cont_Back_ang,600,600);
    clust_continous_angle(i).back_ang = [cont_back_ang.stats];

    cont_Rwrst_dist = [clust_struct_stand_angle(i).clust_rwrst_dist,clust_struct_stand_angle(i).indx];
    cont_Rwrst_dist = sortrows(cont_Rwrst_dist,2); 
    cont_rwrst_dist = rollstat(cont_Rwrst_dist,600,600);
    clust_continous_angle(i).rwrst_dist = [cont_rwrst_dist.stats];

    cont_Lwrst_dist = [clust_struct_stand_angle(i).clust_lwrst_dist,clust_struct_stand_angle(i).indx];
    cont_Lwrst_dist = sortrows(cont_Lwrst_dist,2); 
    cont_lwrst_dist = rollstat(cont_Lwrst_dist,600,600);
    clust_continous_angle(i).lwrst_dist = [cont_lwrst_dist.stats];

end

% For full
full_continous_angle = struct;

clear cont*

cont_SHOxy_ang = [struct_stand_angle.full_SHO_xy_ang,struct_stand_angle.full_indx]; % create matrix of value and index
cont_SHOxy_ang = sortrows(cont_SHOxy_ang,2); % sort according to values
cont_SHO_xy_ang = rollstat(cont_SHOxy_ang,600,600); % calculate rolling window stats (120Hz * 10sec = 600 points)
full_continous_angle.SHO_xy_ang = [cont_SHO_xy_ang.stats]; % save to struct

cont_SHOxz_ang = [struct_stand_angle.full_SHO_xz_ang,struct_stand_angle.full_indx];
cont_SHOxz_ang = sortrows(cont_SHOxz_ang,2); 
cont_SHO_xz_ang = rollstat(cont_SHOxz_ang,600,600);
full_continous_angle.SHO_xz_ang = [cont_SHO_xz_ang.stats];

cont_Head_ang = [struct_stand_angle.full_head_ang,struct_stand_angle.full_indx];
cont_Head_ang = sortrows(cont_Head_ang,2); 
cont_head_ang = rollstat(cont_Head_ang,600,600);
full_continous_angle.head_ang = [cont_head_ang.stats];

cont_Back_ang = [struct_stand_angle.full_back_ang,struct_stand_angle.full_indx];
cont_Back_ang = sortrows(cont_Back_ang,2); 
cont_back_ang = rollstat(cont_Back_ang,600,600);
full_continous_angle.back_ang = [cont_back_ang.stats];

cont_Rwrst_dist = [struct_stand_angle.full_rwrst_dist,struct_stand_angle.full_indx];
cont_Rwrst_dist = sortrows(cont_Rwrst_dist,2); 
cont_rwrst_dist = rollstat(cont_Rwrst_dist,600,600);
full_continous_angle.rwrst_dist = [cont_rwrst_dist.stats];

cont_Lwrst_dist = [struct_stand_angle.full_lwrst_dist,struct_stand_angle.full_indx];
cont_Lwrst_dist = sortrows(cont_Lwrst_dist,2); 
cont_lwrst_dist = rollstat(cont_Lwrst_dist,600,600);
full_continous_angle.lwrst_dist = [cont_lwrst_dist.stats];

disp("    Continous Parameters Calculated");

%% SAVE .MAT FILE
if init_struct.mat_save
    
    save(strcat(p_num,'_clust_stand_angle.mat'),'clust_struct_stand_angle');
    save(strcat(p_num,'_stand_angle.mat'),'struct_stand_angle');

    save(strcat(p_num,'_clust_continous_angle.mat'),'clust_continous_angle');
    save(strcat(p_num,'_full_continous_angle.mat'),'full_continous_angle');

    disp('    .MAT Saved');

end

%% PLOT CONTINOUS PARAMETERS

% cluster continous plots
struct_var_name = fieldnames(clust_continous_angle);
figure;
for i=1:length(struct_var_name)
    for k=1:length(clust_continous_angle)
    
        clust_num = strcat("Cluster ",num2str(k));
        len_interval = length([clust_continous_angle(k).(struct_var_name{i}).mean]);
        
        subplot(2,1,1);
        hold on;
%         bar(linspace(1,len_interval,len_interval),[clust_continous_COP(k).(struct_var_name{i}).mean],'grouped',cVec(k),'DisplayName',clust_num);
        plot([clust_continous_angle(k).(struct_var_name{i}).mean],cVec(k),'DisplayName',clust_num);
        title(strcat(struct_var_name(i),' Mean'),'Interpreter','none');
        ylabel(strcat('5 sec',struct_var_name{i},' mean'),'Interpreter','none'),xlabel('Intervals');
        legend('-DynamicLegend');
        legend('show');
        hold off;

        subplot(2,1,2);
        hold on;
%         bar(linspace(1,len_interval,len_interval),[clust_continous_COP(k).(struct_var_name{i}).sd],'grouped',cVec(k),'DisplayName',clust_num);
        plot([clust_continous_angle(k).(struct_var_name{i}).sd],cVec(k),'DisplayName',clust_num);
        title(strcat(struct_var_name(i),' Std'),'Interpreter','none');
        ylabel(strcat('5 sec',struct_var_name(i),' std'),'Interpreter','none'),xlabel('Intervals');
        legend('-DynamicLegend');
        legend('show');
        hold off;
        
    end
    if init_struct.plot_save % Save if asked
        saveas(gcf,strcat(p_num,'_clust_contionus_',struct_var_name{i},'.jpg'));
    end
    if ~init_struct.plot_view % Keep graph if asked
        close;
    end
    
end

% full continous plots
struct_var_name = fieldnames(full_continous_angle);
figure;
for i=1:length(struct_var_name)
    
    subplot(2,1,1);
    hold on;
    bar([full_continous_angle.(struct_var_name{i}).mean]);
    title(strcat(struct_var_name(i),' Mean'),'Interpreter','none');
    ylabel(strcat('5 sec',struct_var_name{i},' mean'),'Interpreter','none'),xlabel('Intervals');
    hold off;

    subplot(2,1,2);
    hold on;
    bar([full_continous_angle.(struct_var_name{i}).sd]);
    title(strcat(struct_var_name(i),' Std'),'Interpreter','none');
    ylabel(strcat('5 sec',struct_var_name(i),' std'),'Interpreter','none'),xlabel('Intervals');
    hold off;
    
    if init_struct.plot_save % Save if asked
        saveas(gcf,strcat(p_num,'_full_contionus_',struct_var_name{i},'.jpg'));
    end
    if ~init_struct.plot_view % Keep graph if asked
        close;
    end
    
end

disp("    Continous Plots Saved");

%% CHANGE BACK TO DIRECTORY

cd('..')

end