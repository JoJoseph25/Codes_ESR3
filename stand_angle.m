function [lwrst_indx_clust,rwrst_indx_clust,struct_stand_angle] = stand_angle(init_struct,p_num,Mrk_Data,ind_start,ind_end)
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

mrkr_RWRSTz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RWRSTz)/1; % m
mrkr_LWRSTz_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LWRSTz)/1; % m

mrkr_RWRSTy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RWRSTy)/1; % m
mrkr_LWRSTy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LWRSTy)/1; % m


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

stand_LSHO_xz_ang = atan2( stand_LSHO_xz(:,2), stand_LSHO_xz(:,1)).*180./pi; % LSHO C7 ANG
stand_RSHO_xz_ang = atan2( stand_RSHO_xz(:,2), stand_RSHO_xz(:,1)).*180./pi; % RSHO C7 ANG

stand_SHO_xy_ang = stand_LSHO_xy_ang-stand_RSHO_xy_ang; % SHO XY ANG
stand_SHO_xz_ang = stand_LSHO_xz_ang-stand_RSHO_xz_ang; % SHO XZ ANG

stand_head_ang = atan2(stand_head_v2(:,2),stand_head_v2(:,1))-atan2(stand_head_v1(:,2),stand_head_v1(:,1)); % HEAD ANG

stand_back_ang = atan2(stand_back_v2(:,2),stand_back_v2(:,1))-atan2(stand_back_v1(:,2),stand_back_v1(:,1)); %  BACK ANG

disp('    Calculations Done');

%% CLUSTERING
% This is to cluster when you raise hand and keep portion of stable standing

rwrst_clust_mat = horzcat(stand_RWRSTy,stand_RWRSTz); % right wrist matrix
lwrst_clust_mat = horzcat(stand_LWRSTy,stand_LWRSTz); % left wrist matrix

rwrst_clust_mat = transpose(rwrst_clust_mat); % right wrist matrix transpose
lwrst_clust_mat = transpose(lwrst_clust_mat); % left wrist matrix transpose

min_dist=25; % min distance for cluster

% Function of Mean Shift Cluster to form initial clusters on entire standing Right Wrist data
[rwrst_clustCent,rwrst_point2cluster,rwrst_clustMembsCell] = MeanShiftCluster(rwrst_clust_mat,min_dist);
% Function of Mean Shift Cluster to form initial clusters on entire standing Left Wrist data
[lwrst_clustCent,lwrst_point2cluster,lwrst_clustMembsCell] = MeanShiftCluster(lwrst_clust_mat,min_dist);
rwrst_numClust = length(rwrst_clustMembsCell); % number of clusters left wrist
lwrst_numClust = length(lwrst_clustMembsCell); % number of clusters right wrist

%Right Wrist Cleaning
for k=1:rwrst_numClust
    if length(rwrst_clustMembsCell{k})>120
        rwrst_indx_temp = rwrst_clustMembsCell{k}; % current right wrist cluster indexes
        rwrst_dist_temp = rwrst_clust_mat(:,rwrst_clustMembsCell{k}); % current right wrist cluster data points 
        rwrst_dist_temp = transpose(rwrst_dist_temp); % transpose the current right wrist cluster data points

        % create distance vector to be used to clean data or remove transitions
        rwrst_dist_mat = pdist(rwrst_dist_temp,'mahalanobis'); % Malabonis distance vector 
        clear rwrst_dist_temp % remove temp data point matrix as it is not needed
        rwrst_dist_mat = squareform(rwrst_dist_mat); % symmetric rwrist distance matrix

        Eps = 0.72; % min distance 
        MinPts = 120; % min data points need to call as cluster (1 sec=120*1 points)

        % Function of Density based clustering to remove the transitions 
        rwrst_Clust = DBSCAN(rwrst_dist_mat,Eps,MinPts);
        rwrst_indx_clust{k,1} = rwrst_indx_temp(rwrst_Clust~=0); % indexes of clean clusters    
    end
end

%Left Wrist Cleaning
for k=1:lwrst_numClust
    if length(lwrst_clustMembsCell{k})>120
        lwrst_indx_temp = lwrst_clustMembsCell{k}; % current left wrist cluster indexes
        lwrst_dist_temp = lwrst_clust_mat(:,lwrst_clustMembsCell{k}); % current left wrist cluster data points 
        lwrst_dist_temp = transpose(lwrst_dist_temp); % transpose the current left wrist cluster data points

        % create distance vector to be used to clean data or remove transitions
        lwrst_dist_mat = pdist(lwrst_dist_temp,'mahalanobis'); % Malabonis distance vector
        clear lwrst_dist_temp % remove temp data point matrix as it is not needed
        lwrst_dist_mat = squareform(lwrst_dist_mat); % symmetric rwrist distance matrix

        Eps = 0.72; % min distance 
        MinPts = 120; % min data points need to call as cluster (1 sec=120*1 points)

        % Function of Density based clustering to remove the transitions 
        lwrst_Clust = DBSCAN(lwrst_dist_mat,Eps,MinPts);
        lwrst_indx_clust{k,1} = lwrst_indx_temp(lwrst_Clust~=0); % indexes of clean clusters
    end
end

Lwrst_indx = horzcat(lwrst_indx_clust{:}); % Left wrist clean indexes
Rwrst_indx = horzcat(rwrst_indx_clust{:}); % Right wrist clean indexes

wrst_indx = intersect(Lwrst_indx,Rwrst_indx); % Clean Wrist combined indexes

disp('    Clustering');

%% PLOTS

% Unclean
figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
subplot(2,1,1);
hold on;
plot(stand_time,stand_RWRSTz,[cVec(k) '.']);
title("RIGHT WRIST ")
ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;

subplot(2,1,2);
hold on;
plot(stand_time,stand_LWRSTz,[cVec(k) '.']);
title("LEFT WRIST")
ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_wrist_elevation.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Clean
figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
subplot(2,1,1);
hold on;
plot(stand_time(wrst_indx),stand_RWRSTz(wrst_indx),[cVec(k) '.']);
title("RIGHT WRIST Clean")
ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;

subplot(2,1,2);
hold on;
plot(stand_time(wrst_indx),stand_LWRSTz(wrst_indx),[cVec(k) '.']);
title("LEFT WRIST Clean")
ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'clean_wrist_elevation.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Wrist Plots');

%% SHOULDER ANGLES

figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
subplot(2,1,1);
hold on;
plot(stand_time(wrst_indx),stand_SHO_xy_ang(wrst_indx),[cVec(k) '.']);
title("Shoulder XY Angle")
ylabel('Shoulder XY'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;

subplot(2,1,2);
hold on;
plot(stand_time(wrst_indx),stand_SHO_xz_ang(wrst_indx),[cVec(k) '.']);
title("Shoulder XZ Angle")
ylabel('Shoulder XZ'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'shoulder_ang.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

%% NECK & BACK ANGLES

figure;
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
subplot(2,1,1);
hold on;
plot(stand_time(wrst_indx),stand_head_ang(wrst_indx),[cVec(k) '.']);
title("Neck Angle")
ylabel('Neck'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;

subplot(2,1,2);
hold on;
plot(stand_time(wrst_indx),stand_back_ang(wrst_indx),[cVec(k) '.']);
title("Back Angle")
ylabel('Back'),xlabel(strcat('Time (',stand_range,'secs)'));
hold off;
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'back_and_neck_ang.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end


%% SAVE .MAT FILE

struct_stand_angle=struct;

disp('.MAT Saved');
end