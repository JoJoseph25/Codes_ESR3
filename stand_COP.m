function [struct_stand_COP] = stand_COP(init_struct,p_num,Mrk_Data,FP_data,ind_start,ind_end)
% function [struct_stand_COP] = stand_COP(init_struct,p_num,Mrk_Data,FP_data,ind_start,ind_end)
% 
%  This function calculates the COP parameters while saving
%  struct_stand_COP and plots all the graph and .mat file
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Condition Number (string)
%           Mrk_Data - Marker data (table)
%           FP_data - Force plate data (table)
%           ind_start - Start Index to slice data (index)
%           ind_end - End Index to slice data (index)
%     
%   OUTPUT: struct_stand_COP - COP parameters (struct)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% EXCEPTIONS

exceptions={'1008_dy8','1010_dy6','1018_dy7','1018_dy8'}; % List of exceptions that need different 
%                                     min distance than normal data for
%                                     clustering      

disp("COP: ");
%% COP DATA

COPx = FP_data.COP_X(:); % COP X
COPx(isnan(COPx))= 0; % Fill the NaN Values with zero
COPy = FP_data.COP_Y(:); % COP Y
COPy(isnan(COPy))= 0; % Fill the NaN Values with zero

% Filter the data 
COPx_Filter = filtfilt(init_struct.b_fp,init_struct.a_fp,COPx)/1; 
COPy_Filter = filtfilt(init_struct.b_fp,init_struct.a_fp,COPy)/1; 

% COP XY - for the plot
COPxy_Filter = (((COPx_Filter.^2)+(COPy_Filter.^2)).^0.5); 

% COP Standing 
% from 30 seconds - till 5 sec before the subject leaves the FP
stand_range=num2str((ind_end-ind_start)/init_struct.Fs); % Time range 

stand_time = FP_data.TIME(ind_start:ind_end,:); % Timr vector

stand_COP_X = COPx_Filter(ind_start:ind_end); % COP X vector
stand_COP_Y = COPy_Filter(ind_start:ind_end); % COP Y vector
stand_COP_XY= COPxy_Filter(ind_start:ind_end); % COP XY vector

clear COP* % clear unnecessary variables 

%% INITIAL PLOTS

% Plot for COP X vs COP Y for the entire duration
figure;
plot(FP_data.COP_X,FP_data.COP_Y);
title("Full COP X & Y")
xlabel('COP X [mm] (medial-lateral)'),ylabel('COP Y [mm] (anterior-posterior)');
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_fullCOP.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Plot for COP X vs COP Y while struct_stand_COP
figure;
plot(stand_COP_X,stand_COP_Y);
title(strcat('COP X & Y 30:-5 sec before First Step (',stand_range,'secs)'));
xlabel('COP X [mm] (medial-lateral)'),ylabel('COP Y [mm] (anterior-posterior)');
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_standCOP.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Plot for COP X/COP Y vs time while struct_stand_COP
figure;
subplot(2,1,1);
plot(stand_time,stand_COP_X);
title("COP X 30:-5 sec before First Step")
ylabel('COP X [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
subplot(2,1,2);
plot(stand_time,stand_COP_Y);
title("COP X 30:-5 sec before First Step")
ylabel('COP Y [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_COP_time.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Initial Plot');

%% COP CLUSTERING
% This is to cluster COP X-Y data to see wether there was a shift while standing

clust_mat=horzcat(stand_COP_X,stand_COP_Y); % cluster matrix
clust_mat=transpose(clust_mat); % cluster transpose

if ismember(p_num,exceptions) % Check if current subject + trial in exception
    min_dist=15;    % Special Cases minimum distance
else
    min_dist=25;    % Normal minimum distance
end    

% Function of Mean Shift Cluster to form initial clusters on entire standing COP data
[clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(clust_mat,min_dist);

numClust = length(clustMembsCell); % number of clusters

% Clean clusters and remove transitions
for k=1:numClust
    indx_temp=clustMembsCell{k}; % current cluster indexes
    dist_temp=clust_mat(:,clustMembsCell{k}); % current cluster data points (COP X-Y)
    dist_temp=transpose(dist_temp); % transpose the current cluster data points
    
    % create distance vector to be used to clean data or remove transitions
    dist_mat=pdist(dist_temp,'mahalanobis'); % Malabonis distance vector  
    clear dist_temp % remove temp data point matrix as it is not needed
    dist_mat=squareform(dist_mat); % symmetric distance matrix

    Eps=0.72; % min distance 
    MinPts=240; % min data points need to call as cluster (2 sec=120*2 points)
    
    % Function of Density based clustering to remove the transitions 
    Clust = DBSCAN(dist_mat,Eps,MinPts);
    indx_clust{k,1}=indx_temp(Clust~=0); % indexes of clean clusters
end   

disp('    Clustering');

%% CLUSTER PLOTS

% Unclean cluster plot COP x vs COP y
figure;
hold on
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
for k = 1:min(numClust,length(cVec))
    cop_num=strcat("COP",num2str(k));
    cent_num=strcat("Center",num2str(k));
    myMembers = clustMembsCell{k};
    myClustCen = clustCent(:,k);
    plot(clust_mat(1,myMembers),clust_mat(2,myMembers),[cVec(k) '.'],'DisplayName',cop_num);
    plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k),'MarkerSize',10,'DisplayName',cent_num);
    legend('-DynamicLegend');
    legend('show');
end
hold off;
title("COP X & Y")
xlabel('COP X [mm] (medial-lateral)'),ylabel('COP Y [mm] (anterior-posterior)');
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_COP_clust.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end


% Unclean cluster plot COP x/y vs Time
figure;
for k = 1:min(numClust,length(cVec))
    
    subplot(2,1,1);
    hold on;
    plot(stand_time(clustMembsCell{k}),stand_COP_X(clustMembsCell{k}),[cVec(k) '.']);
    title("COP X 30:-5 sec before First Step")
    ylabel('COP X [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off; 

    subplot(2,1,2);
    hold on;
    plot(stand_time(clustMembsCell{k}),stand_COP_Y(clustMembsCell{k}),[cVec(k) '.']);
    title("COP Y 30:-5 sec before First Step")
    ylabel('COP Y [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
    hold off;
end

if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_COP_clust_time.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Cluster Plots');

%% CLEAN CLUSTER PLOTS

% Cluster plot COP x vs COP y without transitions
figure;
hold on
cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
for k = 1:min(numClust,length(cVec))
    cop_num=strcat("COP",num2str(k));
    cent_num=strcat("Center",num2str(k));
    if length(indx_clust{k})>1
        myMembers = indx_clust{k};
        temp_clust=clust_mat(:,myMembers);
        myClustCen = mean(temp_clust,2);
        plot(clust_mat(1,myMembers),clust_mat(2,myMembers),[cVec(k) '.'],'DisplayName',cop_num);
        plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k),'MarkerSize',10,'DisplayName',cent_num);
        legend('-DynamicLegend');
        legend('show');
        clear temp_clust
    end
end
hold off;
title("COP X & Y")
xlabel('COP X [mm] (medial-lateral)'),ylabel('COP Y [mm] (anterior-posterior)');
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_clean_COP_clust.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

% Cluster plot COP x/y vs Time without Transitions
figure;
for k = 1:min(numClust,length(cVec))
    if length(indx_clust{k})>1
        subplot(2,1,1);
        hold on;
        plot(stand_time(indx_clust{k}),stand_COP_X(indx_clust{k}),[cVec(k) '.']);
        title("COP X 30:-5 sec before First Step")
        ylabel('COP X [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
        hold off;

        subplot(2,1,2);
        hold on;
        plot(stand_time(indx_clust{k}),stand_COP_Y(indx_clust{k}),[cVec(k) '.']);
        title("COP Y 30:-5 sec before First Step")
        ylabel('COP Y [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
        hold off;
        
    end
end
if init_struct.plot_save % Save if asked
    saveas(gcf,strcat(p_num,'_clean_COP_clust_time.jpg'));
end
if ~init_struct.plot_view % Keep graph if asked
    close;
end

disp('    Clean Cluster Plots');

%% NUM OF TRANSITION

if numClust>1
    disp("    More than one cluster");
    trans_indx = [];
    Clust_id = [];

    for i=1:numClust
       clear temp_id;
       trans_indx = [trans_indx (clustMembsCell{i})];
       temp_id = ones(1,length(clustMembsCell{i}));
       Clust_id = [Clust_id (temp_id.*i)];

    end
    
    clust_trans = [trans_indx;Clust_id]';
    clust_trans_sort = sortrows(clust_trans);

    trans_mat = abs(diff(clust_trans_sort));

    transitions = sum(trans_mat(:,2)>0);
    
else
    trans_indx = cell2mat(indx_clust);
    disp("    Only one cluster");
   transitions = 0;
   
end

disp(['    Number of Transitions: ',num2str(transitions)]);

%% %% MARKER DATA

% Right Ankle (RANK)
mrkr_RANKx = Mrk_Data.RANKX(:); % RANKX
mrkr_RANKy = Mrk_Data.RANKY(:); % RANKY

% Left Ankle (LANK)
mrkr_LANKx = Mrk_Data.LANKX(:); % LANKX
mrkr_LANKy = Mrk_Data.LANKY(:); % LANKY

% Right Heel (RHEEL)
mrkr_RHEELx = Mrk_Data.RHEELX(:); % RHEELX
mrkr_RHEELy = Mrk_Data.RHEELY(:); % RHEELY

% Left Heel (LHEEL)
mrkr_LHEELx = Mrk_Data.LHEELX(:); % LHEELX
mrkr_LHEELy = Mrk_Data.LHEELY(:); % LHEELY

% Right Toe (RTOE)
mrkr_RTOEx = Mrk_Data.RTOEX(:); % RTOEX
mrkr_RTOEy = Mrk_Data.RTOEY(:); % RTOEY

% Left Toe (LTOE)
mrkr_LTOEx = Mrk_Data.LTOEX(:); % LTOEX
mrkr_LTOEy = Mrk_Data.LTOEY(:); % LTOEY

% Filter Data
mrkr_RANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKx)/1000; % m
mrkr_RANKy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RANKy)/1000; % m 

mrkr_LANKx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKx)/1000; % m
mrkr_LANKy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LANKy)/1000; % m 

mrkr_RHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELx)/1000; % m
mrkr_RHEELy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RHEELy)/1000; % m 

mrkr_LHEELx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELx)/1000; % m
mrkr_LHEELy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LHEELy)/1000; % m 

mrkr_RTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEx)/1000; % m
mrkr_RTOEy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_RTOEy)/1000; % m 

mrkr_LTOEx_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEx)/1000; % m
mrkr_LTOEy_filter = filtfilt(init_struct.b_mrk,init_struct.a_mrk,mrkr_LTOEy)/1000; % m 

% Right Ankle (RANK)
stand_RANKx = mrkr_RANKx_filter(ind_start:ind_end); % RANKX
stand_RANKy = mrkr_RANKy_filter(ind_start:ind_end); % RANKY

% Left Ankle (LANK)
stand_LANKx = mrkr_LANKx_filter(ind_start:ind_end); % LANKX
stand_LANKy = mrkr_LANKy_filter(ind_start:ind_end); % LANKY

% Right HEEL (RHEEL)
stand_RHEELx = mrkr_RHEELx_filter(ind_start:ind_end); % RHEELX
stand_RHEELy = mrkr_RHEELy_filter(ind_start:ind_end); % RHEELY

% Left HEEL (LHEEL)
stand_LHEELx = mrkr_LHEELx_filter(ind_start:ind_end); % LHEELX
stand_LHEELy = mrkr_LHEELy_filter(ind_start:ind_end); % LHEELY

% Right Toe (RTOE)
stand_RTOEx = mrkr_RTOEx_filter(ind_start:ind_end); % RTOEX
stand_RTOEy = mrkr_RTOEy_filter(ind_start:ind_end); % RTOEY

% Right Toe (LTOE)
stand_LTOEx = mrkr_LTOEx_filter(ind_start:ind_end); % LTOEX
stand_LTOEy = mrkr_LTOEy_filter(ind_start:ind_end); % LTOEY

clear mrkr_* % clear unnecessary variables

disp('    Marker Data Fecthed');

%% STEP WIDTH

stand_step_width = sqrt((stand_RANKy - stand_LANKy).^2); 

disp('    Step Width Calculated');

%% BASE OF SUPPORT

xCordinates =[stand_RHEELx,stand_RANKx,stand_RTOEx,stand_LTOEx,stand_LANKx,stand_LHEELx,stand_RHEELx];
yCordinates =[stand_RHEELy,stand_RANKx,stand_RTOEy,stand_LTOEy,stand_LANKy,stand_LHEELy,stand_RHEELy];

stand_BoS = polyarea(xCordinates',yCordinates');

stand_BoS = stand_BoS';

disp('    Base of Support Calculated');

%% COP PARAMETERS 

clust_struct_stand_COP = struct;
struct_stand_COP = struct;

if numClust>1 % More than one cluster
    for k=1:numClust
        if length(indx_clust{k})>1     
            %% Cluster COP Parameters
            
            % Cluster info
            clust_struct_stand_COP(k).indx = indx_clust{k}'; % index
            clust_struct_stand_COP(k).data_pnts = length(indx_clust{k}); % data points
            clust_struct_stand_COP(k).time = length(indx_clust{k})/init_struct.Fs; % time
            
            % COP data
            clust_struct_stand_COP(k).clust_COP_X = stand_COP_X(clust_struct_stand_COP(k).indx)/1000; % COP X of specific cluster in meters
            clust_struct_stand_COP(k).clust_COP_Y = stand_COP_Y(clust_struct_stand_COP(k).indx)/1000; % COP Y of specific cluster in meters
            clust_struct_stand_COP(k).clust_COP_XY = stand_COP_XY(clust_struct_stand_COP(k).indx)/1000; % COP XY of specific cluster
            clust_struct_stand_COP(k).clust_COP_XY_ratio = (clust_struct_stand_COP(k).clust_COP_X)./(clust_struct_stand_COP(k).clust_COP_Y); % COP XY ratio of specific cluster 
            
            % Mean
            clust_struct_stand_COP(k).clust_mean_COP_X = mean(clust_struct_stand_COP(k).clust_COP_X);
            clust_struct_stand_COP(k).clust_mean_COP_Y = mean(clust_struct_stand_COP(k).clust_COP_Y);
            clust_struct_stand_COP(k).clust_mean_COP_XY = mean(clust_struct_stand_COP(k).clust_COP_XY);
            clust_struct_stand_COP(k).clust_mean_COP_XY_ratio = mean(clust_struct_stand_COP(k).clust_COP_XY_ratio);
            
            % Standard-deviation
            clust_struct_stand_COP(k).clust_std_COP_X = std(clust_struct_stand_COP(k).clust_COP_X);
            clust_struct_stand_COP(k).clust_std_COP_Y = std(clust_struct_stand_COP(k).clust_COP_Y);
            clust_struct_stand_COP(k).clust_std_COP_XY = std(clust_struct_stand_COP(k).clust_COP_XY);
            clust_struct_stand_COP(k).clust_std_COP_XY_ratio = std(clust_struct_stand_COP(k).clust_COP_XY_ratio);
            
            % Consequetive element velocity to reduce spikes
            clear temp_diff temp_vel* temp_indx;
            temp_diff = diff(clust_struct_stand_COP(k).indx);
            j=1;
            for i=1:length(temp_diff)
                if temp_diff(i) == 1
                    temp_vel_COPx(j) = (stand_COP_X(clust_struct_stand_COP(k).indx(i+1))-stand_COP_X(clust_struct_stand_COP(k).indx(i)))/(init_struct.dt*1000); % m/sec
                    temp_vel_COPy(j) = (stand_COP_Y(clust_struct_stand_COP(k).indx(i+1))-stand_COP_Y(clust_struct_stand_COP(k).indx(i)))/(init_struct.dt*1000); % m/sec
                    temp_indx(j) = clust_struct_stand_COP(k).indx(i+1);
                    j = j+1;
                end    
            end
            
            % Transpose to row vector
            temp_indx = temp_indx';
            temp_vel_COPx = temp_vel_COPx';
            temp_vel_COPy = temp_vel_COPy';
            
            %Index
            clust_struct_stand_COP(k).clust_vel_indx = temp_indx(2:end);
            
            % Velocity
            clust_struct_stand_COP(k).clust_vel_COP_X = temp_vel_COPx(2:end);
            clust_struct_stand_COP(k).clust_vel_COP_Y = temp_vel_COPy(2:end);
            
            % Velocity Mean
            clust_struct_stand_COP(k).clust_vel_mean_COP_X = mean(abs(temp_vel_COPx));
            clust_struct_stand_COP(k).clust_vel_mean_COP_Y = mean(abs(temp_vel_COPy));
            clust_struct_stand_COP(k).clust_vel_mean_COP_XY = mean(abs(((temp_vel_COPy).^2+(temp_vel_COPy).^2).^0.5));
            clust_struct_stand_COP(k).clust_vel_mean_COP_XY_ratio = mean(abs((temp_vel_COPx)./(temp_vel_COPy)));

            % Velocity Standard-deviation
            clust_struct_stand_COP(k).clust_vel_std_COP_X = std(abs(temp_vel_COPx));
            clust_struct_stand_COP(k).clust_vel_std_COP_Y = std(abs(temp_vel_COPy));
            clust_struct_stand_COP(k).clust_vel_std_COP_XY = std(abs(((temp_vel_COPy).^2+(temp_vel_COPy).^2).^0.5));
            clust_struct_stand_COP(k).clust_vel_std_COP_XY_ratio = std(abs((temp_vel_COPx)./(temp_vel_COPy)));           
            
            
            % Step-width
            clust_struct_stand_COP(k).clust_step_width = stand_step_width(clust_struct_stand_COP(k).indx);
            clust_struct_stand_COP(k).clust_mean_step_width = mean(clust_struct_stand_COP(k).clust_step_width);
            clust_struct_stand_COP(k).clust_std_step_width = std(clust_struct_stand_COP(k).clust_step_width);
            
            % Base of Support
            clust_struct_stand_COP(k).clust_BoS = stand_BoS(clust_struct_stand_COP(k).indx);
            clust_struct_stand_COP(k).clust_mean_BoS = mean(clust_struct_stand_COP(k).clust_BoS);
            clust_struct_stand_COP(k).clust_std_BoS = std(clust_struct_stand_COP(k).clust_BoS);
            
        end
    end
% Only one cluster     
else
    % Cluster info
    clust_struct_stand_COP.indx = trans_indx'; % index
    clust_struct_stand_COP.data_pnts = length(trans_indx); % data points
    clust_struct_stand_COP.time = length(trans_indx)/init_struct.Fs; % time

    % COP data
    clust_struct_stand_COP.clust_COP_X = stand_COP_X(clust_struct_stand_COP.indx)/1000; % COP X of specific cluster in meters
    clust_struct_stand_COP.clust_COP_Y = stand_COP_Y(clust_struct_stand_COP.indx)/1000; % COP Y of specific cluster in meters
    clust_struct_stand_COP.clust_COP_XY = stand_COP_XY(clust_struct_stand_COP.indx)/1000; % COP XY of specific cluster
    clust_struct_stand_COP.clust_COP_XY_ratio = (clust_struct_stand_COP.clust_COP_X)./(clust_struct_stand_COP.clust_COP_Y); % COP XY ratio of specific cluster 

    % Mean
    clust_struct_stand_COP.clust_mean_COP_X = mean(clust_struct_stand_COP.clust_COP_X);
    clust_struct_stand_COP.clust_mean_COP_Y = mean(clust_struct_stand_COP.clust_COP_Y);
    clust_struct_stand_COP.clust_mean_COP_XY = mean(clust_struct_stand_COP.clust_COP_XY);
    clust_struct_stand_COP.clust_mean_COP_XY_ratio = mean(clust_struct_stand_COP.clust_COP_XY_ratio);

    % Standard-deviation
    clust_struct_stand_COP.clust_std_COP_X = std(clust_struct_stand_COP.clust_COP_X);
    clust_struct_stand_COP.clust_std_COP_Y = std(clust_struct_stand_COP.clust_COP_Y);
    clust_struct_stand_COP.clust_std_COP_XY = std(clust_struct_stand_COP.clust_COP_XY);
    clust_struct_stand_COP.clust_std_COP_XY_ratio = std(clust_struct_stand_COP.clust_COP_XY_ratio);

    % Consequetive element velocity to reduce spikes
    clear temp_diff temp_vel* temp_indx;
    temp_diff = diff(clust_struct_stand_COP.indx);
    j=1;
    for i=1:length(temp_diff)
        if temp_diff(i) == 1
            temp_vel_COPx(j) = (stand_COP_X(clust_struct_stand_COP.indx(i+1))-stand_COP_X(clust_struct_stand_COP.indx(i)))/(init_struct.dt*1000); % m/sec
            temp_vel_COPy(j) = (stand_COP_Y(clust_struct_stand_COP.indx(i+1))-stand_COP_Y(clust_struct_stand_COP.indx(i)))/(init_struct.dt*1000); % m/sec
            temp_indx(j) = clust_struct_stand_COP.indx(i+1);
            j = j+1;
        end    
    end
    
    % Transpose to row vector
    temp_indx = temp_indx';
    temp_vel_COPx = temp_vel_COPx';
    temp_vel_COPy = temp_vel_COPy';

    %Index
    clust_struct_stand_COP.clust_vel_indx = temp_indx(2:end);

    % Velocity
    clust_struct_stand_COP.clust_vel_COP_X = temp_vel_COPx(2:end);
    clust_struct_stand_COP.clust_vel_COP_Y = temp_vel_COPy(2:end);

    % Velocity Mean
    clust_struct_stand_COP.clust_vel_mean_COP_X = mean(abs(temp_vel_COPx));
    clust_struct_stand_COP.clust_vel_mean_COP_Y = mean(abs(temp_vel_COPy));
    clust_struct_stand_COP.clust_vel_mean_COP_XY = mean(abs(((temp_vel_COPy).^2+(temp_vel_COPy).^2).^0.5));
    clust_struct_stand_COP.clust_vel_mean_COP_XY_ratio = mean(abs((temp_vel_COPx)./(temp_vel_COPy)));

    % Velocity Standard-deviation
    clust_struct_stand_COP.clust_vel_std_COP_X = std(abs(temp_vel_COPx));
    clust_struct_stand_COP.clust_vel_std_COP_Y = std(abs(temp_vel_COPy));
    clust_struct_stand_COP.clust_vel_std_COP_XY = std(abs(((temp_vel_COPy).^2+(temp_vel_COPy).^2).^0.5));
    clust_struct_stand_COP.clust_vel_std_COP_XY_ratio = std(abs((temp_vel_COPx)./(temp_vel_COPy)));           
    
    % Step-width
    clust_struct_stand_COP.clust_step_width = stand_step_width(clust_struct_stand_COP.indx);
    clust_struct_stand_COP.clust_mean_step_width = mean(clust_struct_stand_COP.clust_step_width);
    clust_struct_stand_COP.clust_std_step_width = std(clust_struct_stand_COP.clust_step_width);

    % Base of Support
    clust_struct_stand_COP.clust_BoS = stand_BoS(clust_struct_stand_COP.indx);
    clust_struct_stand_COP.clust_mean_BoS = mean(clust_struct_stand_COP.clust_BoS);
    clust_struct_stand_COP.clust_std_BoS = std(clust_struct_stand_COP.clust_BoS);
    
end

%% Full COP Parameters

struct_stand_COP.trans_indx = trans_indx';
struct_stand_COP.time = length(trans_indx)/init_struct.Fs; 
struct_stand_COP.transitions = transitions;

% COP data
struct_stand_COP.full_COP_X = stand_COP_X(trans_indx)/1000; % COP X in meters
struct_stand_COP.full_COP_Y = stand_COP_Y(trans_indx)/1000; % COP Y in meters
struct_stand_COP.full_COP_XY = stand_COP_XY(trans_indx)/1000; % COP XY 
struct_stand_COP.full_COP_XY_ratio = (struct_stand_COP.full_COP_X)./(struct_stand_COP.full_COP_Y); % COP XY ratio

% Mean
struct_stand_COP.full_mean_COP_X = mean(struct_stand_COP.full_COP_X);
struct_stand_COP.full_mean_COP_Y = mean(struct_stand_COP.full_COP_Y);
struct_stand_COP.full_mean_COP_XY = mean(struct_stand_COP.full_COP_XY);
struct_stand_COP.full_mean_COP_XY_ratio = mean(struct_stand_COP.full_COP_XY_ratio);

% Standard-deviation
struct_stand_COP.full_std_COP_X = std(struct_stand_COP.full_COP_X);
struct_stand_COP.full_std_COP_Y = std(struct_stand_COP.full_COP_Y);
struct_stand_COP.full_std_COP_XY = std(struct_stand_COP.full_COP_XY);
struct_stand_COP.full_std_COP_XY_ratio = std(struct_stand_COP.full_COP_XY_ratio);

% Consequetive element velocity to reduce spikes
clear temp_diff temp_vel* temp_indx;
temp_diff = diff(trans_indx);
j=1;
for i=1:length(temp_diff)
    if temp_diff(i) == 1
        temp_vel_COPx(j) = (stand_COP_X(trans_indx(i+1))-stand_COP_X(trans_indx(i)))/(init_struct.dt*1000); % m
        temp_vel_COPy(j) = (stand_COP_Y(trans_indx(i+1))-stand_COP_Y(trans_indx(i)))/(init_struct.dt*1000); % m
        temp_indx(j) = trans_indx(i+1);
        j = j+1;
    end    
end

% Transpose to row vector
temp_indx = temp_indx';
temp_vel_COPx = temp_vel_COPx';
temp_vel_COPy = temp_vel_COPy';

%Index
struct_stand_COP.full_vel_indx = temp_indx(2:end);

% Velocity
struct_stand_COP.full_vel_COP_X = temp_vel_COPx(2:end);
struct_stand_COP.full_vel_COP_Y = temp_vel_COPy(2:end);

% Velocity Mean 
struct_stand_COP.full_vel_mean_COP_X = mean(abs(temp_vel_COPx));
struct_stand_COP.full_vel_mean_COP_Y = mean(abs(temp_vel_COPy));
struct_stand_COP.full_vel_mean_COP_XY = mean(abs(((temp_vel_COPy).^2+(temp_vel_COPy).^2).^0.5));
struct_stand_COP.full_vel_mean_COP_XY_ratio = mean(abs((temp_vel_COPx)./(temp_vel_COPy)));

% Velocity Standard-deviation 
struct_stand_COP.full_vel_std_COP_X = std(abs(temp_vel_COPx));
struct_stand_COP.full_vel_std_COP_Y = std(abs(temp_vel_COPy));
struct_stand_COP.full_vel_std_COP_XY = std(abs(((temp_vel_COPy).^2+(temp_vel_COPy).^2).^0.5));
struct_stand_COP.full_vel_std_COP_XY_ratio = std(abs((temp_vel_COPx)./(temp_vel_COPy)));           

% Step-width
struct_stand_COP.full_step_width = stand_step_width(trans_indx);
struct_stand_COP.full_mean_step_width = mean(struct_stand_COP.full_step_width);
struct_stand_COP.full_std_step_width = std(struct_stand_COP.full_step_width);

% Base of Support
struct_stand_COP.full_BoS = stand_BoS(trans_indx);
struct_stand_COP.full_mean_BoS = mean(struct_stand_COP.full_BoS);
struct_stand_COP.full_std_BoS = std(struct_stand_COP.full_BoS);

%% RELATIVE PARAMETERS

% parameters adjusted by cluster-size (n*P/n) where n= num of points and
%                                                   P= parameter of interest

% Relative Mean
struct_stand_COP.clust_rel_mean_COP_X = (sum([clust_struct_stand_COP.clust_mean_COP_X].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_mean_COP_Y = (sum([clust_struct_stand_COP.clust_mean_COP_Y].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_mean_COP_XY = (sum([clust_struct_stand_COP.clust_mean_COP_XY].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_mean_COP_XY_ratio = (sum([clust_struct_stand_COP.clust_mean_COP_XY_ratio].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);

% Relative Standard-deviation
struct_stand_COP.clust_rel_std_COP_X = (sum([clust_struct_stand_COP.clust_std_COP_X].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_stsd_COP_Y = (sum([clust_struct_stand_COP.clust_std_COP_Y].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_std_COP_XY = (sum([clust_struct_stand_COP.clust_std_COP_XY].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_std_COP_XY_ratio = (sum([clust_struct_stand_COP.clust_std_COP_XY_ratio].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);

% Relative Velocity Mean
struct_stand_COP.clust_rel_vel_mean_COP_X = (sum([clust_struct_stand_COP.clust_vel_mean_COP_X].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_vel_mean_COP_Y = (sum([clust_struct_stand_COP.clust_vel_mean_COP_Y].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_vel_mean_COP_XY = (sum([clust_struct_stand_COP.clust_vel_mean_COP_XY].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_vel_mean_COP_XY_ratio = (sum([clust_struct_stand_COP.clust_vel_mean_COP_XY_ratio].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);   

% Relative Velocity Standard-deviation
struct_stand_COP.clust_rel_vel_std_COP_X = (sum([clust_struct_stand_COP.clust_vel_std_COP_X].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_vel_std_COP_Y = (sum([clust_struct_stand_COP.clust_vel_std_COP_Y].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_vel_std_COP_XY = (sum([clust_struct_stand_COP.clust_vel_std_COP_XY].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_vel_std_COP_XY_ratio = (sum([clust_struct_stand_COP.clust_vel_std_COP_XY_ratio].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);   

% Step-width
struct_stand_COP.clust_rel_mean_step_width = (sum([clust_struct_stand_COP.clust_mean_step_width].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_std_step_width = (sum([clust_struct_stand_COP.clust_std_step_width].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);

% Base of Support
struct_stand_COP.clust_rel_mean_BoS = (sum([clust_struct_stand_COP.clust_mean_BoS].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);
struct_stand_COP.clust_rel_std_BoS = (sum([clust_struct_stand_COP.clust_std_BoS].*[clust_struct_stand_COP.data_pnts]))/sum([clust_struct_stand_COP.data_pnts]);;

disp("    Parameters Calculated");

%% CONTINOUS PARAMETERS

% For clusters
clust_continous_COP=struct;

for i=1:length(clust_struct_stand_COP)
   clear COP*
   
   COPX = [clust_struct_stand_COP(i).clust_COP_X,clust_struct_stand_COP(i).indx]; % create matrix of value and index
   COPX = sortrows(COPX,2); % sort according to values
   COP_X = rollstat(COPX,360,360); % calculate rolling window stats (120Hz * 3sec = 360 points)
   clust_continous_COP(i).COP_X = [COP_X.stats]; % save to struct
   
   COPY = [clust_struct_stand_COP(i).clust_COP_Y,clust_struct_stand_COP(i).indx];
   COPY = sortrows(COPY,2);
   COP_Y = rollstat(COPY,360,360);
   clust_continous_COP(i).COP_Y = [COP_Y.stats];
   
   COPXY = [clust_struct_stand_COP(i).clust_COP_XY,clust_struct_stand_COP(i).indx];
   COPXY = sortrows(COPXY,2);
   COP_XY = rollstat(COPXY,360,360);
   clust_continous_COP(i).COP_XY = [COP_XY.stats];
   
   COPXY_ratio = [clust_struct_stand_COP(i).clust_COP_XY,clust_struct_stand_COP(i).indx];
   COPXY_ratio = sortrows(COPXY_ratio,2);
   COP_XY_ratio = rollstat(COPXY_ratio,360,360);
   clust_continous_COP(i).COP_XY_ratio = [COP_XY_ratio.stats];
   
   COPX_vel = [clust_struct_stand_COP(i).clust_vel_COP_X,clust_struct_stand_COP(i).clust_vel_indx];
   COPX_vel = sortrows(COPX_vel,2);
   COP_X_vel = rollstat(COPX_vel,360,360);
   clust_continous_COP(i).COP_X_vel = [COP_X_vel.stats];
   
   COPY_vel = [clust_struct_stand_COP(i).clust_vel_COP_Y,clust_struct_stand_COP(i).clust_vel_indx];
   COPY_vel = sortrows(COPY_vel,2);
   COP_Y_vel = rollstat(COPY_vel,360,360);
   clust_continous_COP(i).COP_Y_vel = [COP_Y_vel.stats];
   
   stp_wdt = [clust_struct_stand_COP(i).clust_step_width,clust_struct_stand_COP(i).indx];
   stp_wdt = sortrows(stp_wdt,2);
   step_width = rollstat(stp_wdt,360,360);
   clust_continous_COP(i).step_width = [step_width.stats];
   
   B_o_S = [clust_struct_stand_COP(i).clust_BoS,clust_struct_stand_COP(i).indx];
   B_o_S = sortrows(B_o_S,2);
   BoS = rollstat(B_o_S,360,360);
   clust_continous_COP(i).BoS = [BoS.stats];
  
end    

% For full
full_continous_COP = struct;

clear COP*

COPX = [struct_stand_COP.full_COP_X,struct_stand_COP.trans_indx]; % create matrix of value and index
COPX = sortrows(COPX,2); % sort according to values
COP_X = rollstat(COPX,360,360); % calculate rolling window stats
full_continous_COP.COP_X = [COP_X.stats]; % save to struct

COPY = [struct_stand_COP.full_COP_Y,struct_stand_COP.trans_indx];
COPY = sortrows(COPY,2);
COP_Y = rollstat(COPY,360,360);
full_continous_COP.COP_Y = [COP_Y.stats];

COPX_vel = [struct_stand_COP.full_vel_COP_X,struct_stand_COP.full_vel_indx];
COPX_vel = sortrows(COPX_vel,2);
COP_X_vel = rollstat(COPX_vel,360,360);
full_continous_COP.COP_X_vel = [COP_X_vel.stats];

COPY_vel = [struct_stand_COP.full_vel_COP_Y,struct_stand_COP.full_vel_indx];
COPY_vel = sortrows(COPY_vel,2);
COP_Y_vel = rollstat(COPY_vel,360,360);
full_continous_COP.COP_Y_vel = [COP_Y_vel.stats];   

stp_wdt = [struct_stand_COP.full_step_width,struct_stand_COP.trans_indx];
stp_wdt = sortrows(stp_wdt,2);
step_width = rollstat(stp_wdt,360,360);
full_continous_COP.step_width = [step_width.stats];
   
B_o_S = [struct_stand_COP.full_BoS,struct_stand_COP.trans_indx];
B_o_S = sortrows(B_o_S,2);
BoS = rollstat(B_o_S,360,360);
full_continous_COP.BoS = [BoS.stats];

disp("    Continous Parameters Calculated");

%% MAT FILE

save(strcat(p_num,'_clust_stand_COP.mat'),'clust_struct_stand_COP');
save(strcat(p_num,'_stand_COP.mat'),'struct_stand_COP');

save(strcat(p_num,'_clust_continous_COP.mat'),'clust_continous_COP');
save(strcat(p_num,'_full_continous_COP.mat'),'full_continous_COP');

disp("    COP .MAT File Saved");

%% PLOT CONTINOUS PARAMETERS

% cluster continous plots
struct_var_name = fieldnames(clust_continous_COP);
figure;
for i=1:length(struct_var_name)
    for k=1:length(clust_continous_COP)
    
        clust_num = strcat("Cluster ",num2str(k));
        len_interval = length([clust_continous_COP(k).(struct_var_name{i}).mean]);
        
        subplot(2,1,1);
        hold on;
%         bar(linspace(1,len_interval,len_interval),[clust_continous_COP(k).(struct_var_name{i}).mean],'grouped',cVec(k),'DisplayName',clust_num);
        plot([clust_continous_COP(k).(struct_var_name{i}).mean],cVec(k),'DisplayName',clust_num);
        title(strcat(struct_var_name(i),' Mean'),'Interpreter','none');
        ylabel(strcat('3 sec',struct_var_name{i},' mean'),'Interpreter','none'),xlabel('Intervals');
        legend('-DynamicLegend');
        legend('show');
        hold off;

        subplot(2,1,2);
        hold on;
%         bar(linspace(1,len_interval,len_interval),[clust_continous_COP(k).(struct_var_name{i}).sd],'grouped',cVec(k),'DisplayName',clust_num);
        plot([clust_continous_COP(k).(struct_var_name{i}).sd],cVec(k),'DisplayName',clust_num);
        title(strcat(struct_var_name(i),' Std'),'Interpreter','none');
        ylabel(strcat('3 sec',struct_var_name(i),' std'),'Interpreter','none'),xlabel('Intervals');
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
struct_var_name = fieldnames(full_continous_COP);
figure;
for i=1:length(struct_var_name)
    
    subplot(2,1,1);
    hold on;
    bar([full_continous_COP.(struct_var_name{i}).mean]);
    title(strcat(struct_var_name(i),' Mean'),'Interpreter','none');
    ylabel(strcat('3 sec',struct_var_name{i},' mean'),'Interpreter','none'),xlabel('Intervals');
    hold off;

    subplot(2,1,2);
    hold on;
    bar([full_continous_COP.(struct_var_name{i}).sd]);
    title(strcat(struct_var_name(i),' Std'),'Interpreter','none');
    ylabel(strcat('3 sec',struct_var_name(i),' std'),'Interpreter','none'),xlabel('Intervals');
    hold off;
    
    if init_struct.plot_save % Save if asked
        saveas(gcf,strcat(p_num,'_full_contionus_',struct_var_name{i},'.jpg'));
    end
    if ~init_struct.plot_view % Keep graph if asked
        close;
    end
    
end

disp("    Continous Plots Saved");


end