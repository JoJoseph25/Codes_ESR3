function [struct_stand_COP] = stand_COP(init_struct,p_num,FP_data,ind_start,ind_end)
% function [struct_stand_COP] = stand_COP(init_struct,p_num,FP_data,ind_start,ind_end)
% 
%  This function calculates the COP parameters while saving
%  struct_stand_COP and plots all the graph and .mat file
%
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%           p_num - Subject ID + Experiment Condition Number (string)
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
       trans_indx = [trans_indx (indx_clust{i})];
       temp_id = ones(1,length(indx_clust{i}));
       Clust_id = [Clust_id (temp_id.*i)];

    end
    
    clust_trans = [trans_indx;Clust_id]';
    clust_trans_sort = sortrows(clust_trans);

    trans_mat = diff(clust_trans_sort);

    transitions = sum(trans_mat(:,1)>1 | trans_mat(:,2)>0);
    
else
    trans_indx = cell2mat(indx_clust);
    disp("    Only one cluster");
   transitions = 0;
end

disp(['    Number of Transitions: ',num2str(transitions)]);

%% COP PARAMETERS 

clust_struct_stand_COP=struct;
struct_stand_COP=struct;

if numClust>1 % More than one cluster
    for k=1:numClust
        if length(indx_clust{k})>1     
            %% Cluster COP Parameters
            
            % Cluster info
            clust_struct_stand_COP(k).indx = sort(indx_clust{k}); % index
            clust_struct_stand_COP(k).data_pnts = length(indx_clust{k}); % data points
            clust_struct_stand_COP(k).time = length(indx_clust{k})/init_struct.Fs; % time
            
            % COP data
            clust_struct_stand_COP(k).clust_COP_X = stand_COP_X(clust_struct_stand_COP(k).indx)/1000; % COP X of specific cluster in meters
            clust_struct_stand_COP(k).clust_COP_Y = stand_COP_Y(clust_struct_stand_COP(k).indx)/1000; % COP Y of specific cluster in meters
            clust_struct_stand_COP(k).clust_COP_XY = (((clust_struct_stand_COP(k).clust_COP_X).^2)+((clust_struct_stand_COP(k).clust_COP_X).^2)).^0.5; % COP XY of specific cluster
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
            clear temp_diff temp_vel*;
            temp_diff = diff(clust_struct_stand_COP(k).indx);
            j=1;
            for i=1:length(temp_diff)
                if temp_diff(i) == 1
                    temp_vel_COPx(j) = (stand_COP_X(clust_struct_stand_COP(k).indx(i+1))-stand_COP_X(clust_struct_stand_COP(k).indx(i)))/(init_struct.dt*1000); % m/sec
                    temp_vel_COPy(j) = (stand_COP_Y(clust_struct_stand_COP(k).indx(i+1))-stand_COP_Y(clust_struct_stand_COP(k).indx(i)))/(init_struct.dt*1000); % m/sec
                    j = j+1;
                end    
            end
            
            % Velocity
            clust_struct_stand_COP(k).clust_vel_COP_X = temp_vel_COPx;
            clust_struct_stand_COP(k).clust_vel_COP_Y = temp_vel_COPy;
            
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

        end
    end
% Only one cluster     
else
    % Cluster info
    clust_struct_stand_COP.indx = sort(trans_indx); % index
    clust_struct_stand_COP.data_pnts = length(trans_indx); % data points
    clust_struct_stand_COP.time = length(trans_indx)/init_struct.Fs; % time

    % COP data
    clust_struct_stand_COP.clust_COP_X = stand_COP_X(clust_struct_stand_COP.indx)/1000; % COP X of specific cluster in meters
    clust_struct_stand_COP.clust_COP_Y = stand_COP_Y(clust_struct_stand_COP.indx)/1000; % COP Y of specific cluster in meters
    clust_struct_stand_COP.clust_COP_XY = (((clust_struct_stand_COP.clust_COP_X).^2)+((clust_struct_stand_COP.clust_COP_X).^2)).^0.5; % COP XY of specific cluster
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
    clear temp_diff temp_vel*;
    temp_diff = diff(clust_struct_stand_COP.indx);
    j=1;
    for i=1:length(temp_diff)
        if temp_diff(i) == 1
            temp_vel_COPx(j) = (stand_COP_X(clust_struct_stand_COP.indx(i+1))-stand_COP_X(clust_struct_stand_COP.indx(i)))/(init_struct.dt*1000); % m/sec
            temp_vel_COPy(j) = (stand_COP_Y(clust_struct_stand_COP.indx(i+1))-stand_COP_Y(clust_struct_stand_COP.indx(i)))/(init_struct.dt*1000); % m/sec
            j = j+1;
        end    
    end

    % Velocity
    clust_struct_stand_COP.clust_vel_COP_X = temp_vel_COPx;
    clust_struct_stand_COP.clust_vel_COP_Y = temp_vel_COPy;

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
    
end

%% Full COP Parameters

struct_stand_COP.trans_indx = sort(trans_indx);
struct_stand_COP.time = length(trans_indx)/init_struct.Fs; 
struct_stand_COP.transitions = transitions;

% COP data
struct_stand_COP.full_COP_X = stand_COP_X(trans_indx)/1000; % COP X of specific cluster in meters
struct_stand_COP.full_COP_Y = stand_COP_Y(trans_indx)/1000; % COP Y of specific cluster in meters
struct_stand_COP.full_COP_XY = (((struct_stand_COP.full_COP_X).^2)+((struct_stand_COP.full_COP_X).^2)).^0.5; % COP XY of specific cluster
struct_stand_COP.full_COP_XY_ratio = (struct_stand_COP.full_COP_X)./(struct_stand_COP.full_COP_Y); % COP XY ratio of specific cluster

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
clear temp_diff temp_vel*;
temp_diff = diff(trans_indx);
j=1;
for i=1:length(temp_diff)
    if temp_diff(i) == 1
        temp_vel_COPx(j) = (stand_COP_X(trans_indx(i+1))-stand_COP_X(trans_indx(i)))/(init_struct.dt*1000); % m
        temp_vel_COPy(j) = (stand_COP_Y(trans_indx(i+1))-stand_COP_Y(trans_indx(i)))/(init_struct.dt*1000); % m
        j = j+1;
    end    
end

% Velocity
struct_stand_COP.full_vel_COP_X = temp_vel_COPx;
struct_stand_COP.full_vel_COP_Y = temp_vel_COPy;

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

disp("    Parameters Calculated");

%% MAT FILE

save(strcat(p_num,'_clust_stand_COP.mat'),'clust_struct_stand_COP');
save(strcat(p_num,'_stand_COP.mat'),'struct_stand_COP');

disp("    .MAT File Saved");

end