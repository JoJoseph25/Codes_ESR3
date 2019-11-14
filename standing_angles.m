clear all; clc; close all;

main_workdir = 'G:\Socrates\Codes'; % Lab Computer  

cd(main_workdir);

exp_workdir='G:\Socrates\Exp_Data\Henrique_data\Exp_data';

subject_num={'1017'};
exp_id={'dy7'};


%% PATH CHNG
for sub=1:length(subject_num)
    for exp=1:length(exp_id)
        exp
        cd(main_workdir);

        x_num=subject_num{sub};

        sub_file=strcat(exp_workdir,'\',x_num);

        addpath(main_workdir)
        addpath(sub_file)

        cd(sub_file)

        p_num=strcat(x_num,'_',exp_id{exp});

%% Lowpass Butterworth filter 

        % Force Plate Data -> cutoff frequency 20 Hz
        Fc = 20; % cuttoff frequency
        Fs = 120; % sampling rate
        order = 2; % order
        [b_fp,a_fp] = butter(order,2*Fc/Fs);

        % Markers Data -> cutoff frequency 10 Hz
        Fc_mrk = 10; % cuttoff frequency
        [b_mrk,a_mrk] = butter(order,2*Fc_mrk/Fs);

        % File Names
        filename_mrk = strcat(p_num,'.tsv'); % Force Plate
        filename_FP = strcat(p_num,'_f_1.tsv'); % First Step + Upper Body 

%% IMPORTING
        % Import File Force Plate 
        opts_FP = detectImportOptions(filename_FP,'FileType','text');
        FP_data = readtable(filename_FP,opts_FP);

        % Import File First Step
        opts_mrk = detectImportOptions(filename_mrk,'FileType','text');
        opts_mrk.VariableNamesLine=11;
        opts_mrk.DataLine=12;
        Mrk_Data = readtable(filename_mrk,opts_mrk,'ReadVariableNames',true);
        t_time_frst_stpep = Mrk_Data.Time(:);

        rmpath(sub_file)

%% Quiet Standing Markers

        % Time
        f30_sec=30*120; % Remove first 30-sec 
        dt = 1/Fs; % Fs=120Hz
        t = FP_data.TIME(:); % Time

        % Time point when the subject moves from Force Plate

        % Value of Force Z around zero (first value Force_Z>-5 + RHEEL\LHEEL>2)
        Force_Z = FP_data.Force_Z(:); % Force Z
        RHEEL_X = Mrk_Data.RHEELX(:); % RHEEL X
        LHEEL_X = Mrk_Data.LHEELX(:); % LHEEL X
        
        indx_frst_stp = find((Force_Z>-5 & (RHEEL_X >50 | LHEEL_X > 50)),1,'first'); % Indicator First Step Initiation
        time_frst_stp = t(indx_frst_stp); % Time First Step Initiation
        
        % For COP X, COP X Speed Difference and First Step
        indx_bfr_frst_stp = 3600; % for the COP max. min. changes; 30 seconds before that 30*120 = 3600
        indx_frst_stp_strt = indx_frst_stp-indx_bfr_frst_stp; 
        indx_frst_stp_end = indx_frst_stp-600; % -600 (5 sec) in order to avoid the max. or min. %close to the point where subject leaves the FP

        % Additional Indicators 
        min_time = min(t); % Start Time
        max_time = max(t); % End Time
        t_len = max_time - min_time; % Time Length 

        ssfs = time_frst_stp-min_time; % Seconds from Start to First Step 
        sfse = max_time-time_frst_stp; % Seconds from First Step till End 

        % Set the Time Range
        r1 = 10; % Start Time for Natural (10 seconds) r like range
        r2 = t(indx_frst_stp); % Time range end point [sec]

        t_rng = r2 - r1; % Everything will be calculated for time range

        % Define the start-end rows for the time range

        % Define the Start
        ind_start = f30_sec; % data after first 30 seconds

        % Define the End
        ind_end=indx_frst_stp_end; %5 sec before first step initiation

        % Markers
        TPHDx = Mrk_Data.TPHDX(:); % TPHDX
        TPHDy = Mrk_Data.TPHDY(:); % TPHDY
        TPHDz = Mrk_Data.TPHDZ(:); % TPHDZ
        
        C7x = Mrk_Data.C7X(:); % C7X
        C7y = Mrk_Data.C7Y(:); % C7Y
        C7z = Mrk_Data.C7Z(:); % C7Z
        
        RSHOx = Mrk_Data.RSHOX(:); % RSHOX
        RSHOy = Mrk_Data.RSHOY(:); % RSHOY
        RSHOz = Mrk_Data.RSHOZ(:); % RSHOZ
        
        LSHOx = Mrk_Data.LSHOX(:); % LSHOX
        LSHOy = Mrk_Data.LSHOY(:); % LSHOY
        LSHOz = Mrk_Data.LSHOZ(:); % LSHOZ
        
        IPx = Mrk_Data.IPX(:); % IPX
        IPy = Mrk_Data.IPY(:); % IPY
        IPz = Mrk_Data.IPZ(:); % IPZ
        
        RWRSTz = Mrk_Data.RWRSTZ(:); % RWRSTZ
        LWRSTz = Mrk_Data.LWRSTZ(:); % LWRSTZ
        
        RWRSTy = Mrk_Data.RWRSTY(:); % RWRSTX
        LWRSTy = Mrk_Data.LWRSTY(:); % LWRSTX
        
        num_rec = length(TPHDx); % Number of records 

        % Use the filter
        TPHDx_filter = filtfilt(b_mrk,a_mrk,TPHDx)/1000; % m
        TPHDy_filter = filtfilt(b_mrk,a_mrk,TPHDy)/1000; % m % TPHDY
        TPHDz_filter = filtfilt(b_mrk,a_mrk,TPHDz)/1000; % m
        
        C7x_filter = filtfilt(b_mrk,a_mrk,C7x)/1000; % m
        C7y_filter = filtfilt(b_mrk,a_mrk,C7y)/1000; % m
        C7z_filter = filtfilt(b_mrk,a_mrk,C7z)/1000; % m
        
        RSHOx_filter = filtfilt(b_mrk,a_mrk,RSHOx)/1000; % m
        RSHOy_filter = filtfilt(b_mrk,a_mrk,RSHOy)/1000; % m
        RSHOz_filter = filtfilt(b_mrk,a_mrk,RSHOz)/1000; % m
        
        LSHOx_filter = filtfilt(b_mrk,a_mrk,LSHOx)/1000; % m
        LSHOy_filter = filtfilt(b_mrk,a_mrk,LSHOy)/1000; % m
        LSHOz_filter = filtfilt(b_mrk,a_mrk,LSHOz)/1000; % m
        
        IPx_filter = filtfilt(b_mrk,a_mrk,IPx)/1000; % m
        IPy_filter = filtfilt(b_mrk,a_mrk,IPy)/1000; % m
        IPz_filter = filtfilt(b_mrk,a_mrk,IPz)/1000; % m
        
        RWRSTz_filter = filtfilt(b_mrk,a_mrk,RWRSTz)/1; % m
        LWRSTz_filter = filtfilt(b_mrk,a_mrk,LWRSTz)/1; % m
        
        RWRSTy_filter = filtfilt(b_mrk,a_mrk,RWRSTy)/1; % m
        LWRSTy_filter = filtfilt(b_mrk,a_mrk,LWRSTy)/1; % m
        
        % COP Quiet Standing 
        % from 30 seconds - till 5 sec before the subject leaves the FP
        stand_range=num2str((ind_end-ind_start)/Fs);

        stand_time = Mrk_Data.Time(ind_start:ind_end,:);

        stand_TPHDx = TPHDx_filter(ind_start:ind_end);
        stand_TPHDy = TPHDy_filter(ind_start:ind_end);
        stand_TPHDz = TPHDz_filter(ind_start:ind_end);
        
        stand_C7x = C7x_filter(ind_start:ind_end);
        stand_C7y = C7y_filter(ind_start:ind_end);
        stand_C7z = C7z_filter(ind_start:ind_end);
        
        stand_RSHOx = RSHOx_filter(ind_start:ind_end);
        stand_RSHOy = RSHOy_filter(ind_start:ind_end);
        stand_RSHOz = RSHOz_filter(ind_start:ind_end);
        
        stand_LSHOx = LSHOx_filter(ind_start:ind_end);
        stand_LSHOy = LSHOy_filter(ind_start:ind_end);
        stand_LSHOz = LSHOz_filter(ind_start:ind_end);
        
        stand_IPx = IPx_filter(ind_start:ind_end);
        stand_IPy = IPy_filter(ind_start:ind_end);
        stand_IPz = IPz_filter(ind_start:ind_end);
        
        stand_RWRSTz = RWRSTz_filter(ind_start:ind_end);
        stand_LWRSTz = LWRSTz_filter(ind_start:ind_end);
        
        stand_RWRSTy = RWRSTy_filter(ind_start:ind_end);
        stand_LWRSTy = LWRSTy_filter(ind_start:ind_end);
        
        clearvars -except exp sub_num main_workdir exp_workdir subject_num exp_id...
            stand_range stand_time stand_TPHDx stand_TPHDy stand_TPHDz ...
            stand_C7x stand_C7y stand_C7z stand_RSHOx stand_RSHOy stand_RSHOz ...
            stand_LSHOx stand_LSHOy stand_LSHOz stand_IPx stand_IPy stand_IPz ...
            stand_RWRSTz stand_LWRSTz stand_RWRSTy stand_LWRSTy
        
%% Angle Calculation
        stand_RSHO_xy=[stand_RSHOx-stand_C7x stand_RSHOy-stand_C7y];
        stand_LSHO_xy=[stand_LSHOx-stand_C7x stand_LSHOy-stand_C7y];
        
        stand_RSHO_xz=[stand_RSHOx-stand_C7x stand_RSHOz-stand_C7z];
        stand_LSHO_xz=[stand_LSHOx-stand_C7x stand_LSHOz-stand_C7z];
        
        stand_head_v1=[stand_TPHDx-stand_C7x stand_TPHDz-stand_C7z];
        stand_head_v2=[(stand_TPHDx-stand_C7x)*0 stand_TPHDz-stand_C7z];
        
        stand_back_v1=[stand_IPx-stand_C7x stand_IPz-stand_C7z];
        stand_back_v2=[(stand_IPx-stand_C7x)*0 stand_IPz-stand_C7z];
        
        % Upper Body Angles
        stand_LSHO_xy_ang = atan2( stand_LSHO_xy(:,2), stand_LSHO_xy(:,1)).*180./pi; % Lsho C7
        stand_RSHO_xy_ang = atan2( stand_RSHO_xy(:,2), stand_RSHO_xy(:,1)).*180./pi; % Rsho C7
        
        stand_LSHO_xz_ang = atan2( stand_LSHO_xz(:,2), stand_LSHO_xz(:,1)).*180./pi; % Lsho C7
        stand_RSHO_xz_ang = atan2( stand_RSHO_xz(:,2), stand_RSHO_xz(:,1)).*180./pi; % Rsho C7
        
        stand_SHO_xy_ang = stand_LSHO_xy_ang-stand_RSHO_xy_ang;
        stand_SHO_xz_ang = stand_LSHO_xz_ang-stand_RSHO_xz_ang;
        
        stand_head_ang = atan2(stand_head_v2(:,2),stand_head_v2(:,1))-atan2(stand_head_v1(:,2),stand_head_v1(:,1));
        
        stand_back_ang = atan2(stand_back_v2(:,2),stand_back_v2(:,1))-atan2(stand_back_v1(:,2),stand_back_v1(:,1));
        
 %% CLUSTERING
        
        rwrst_clust_mat=horzcat(stand_RWRSTy,stand_RWRSTz);
        lwrst_clust_mat=horzcat(stand_LWRSTy,stand_LWRSTz);
        
        rwrst_clust_mat=transpose(rwrst_clust_mat);
        lwrst_clust_mat=transpose(lwrst_clust_mat);

        min_dist=25;

        [rwrst_clustCent,rwrst_point2cluster,rwrst_clustMembsCell] = MeanShiftCluster(rwrst_clust_mat,min_dist);
        [lwrst_clustCent,lwrst_point2cluster,lwrst_clustMembsCell] = MeanShiftCluster(lwrst_clust_mat,min_dist);
        
        rwrst_numClust = length(rwrst_clustMembsCell);
        lwrst_numClust = length(lwrst_clustMembsCell);
        
        %Right Wrist Cleaning
        for k=1:rwrst_numClust
            rwrst_indx_temp=rwrst_clustMembsCell{k};
            rwrst_dist_temp=rwrst_clust_mat(:,rwrst_clustMembsCell{k});
            rwrst_dist_temp=transpose(rwrst_dist_temp);
            rwrst_dist_mat=pdist(rwrst_dist_temp,'mahalanobis');
            rwrst_dist_mat=squareform(rwrst_dist_mat);
            
            Eps=0.72;
            MinPts=120;

            rwrst_Clust = DBSCAN(rwrst_dist_mat,Eps,MinPts);
            rwrst_indx_clust{k,1}=rwrst_indx_temp(rwrst_Clust~=0);
        end
        
        %Left Wrist Cleaning
        for k=1:lwrst_numClust
            lwrst_indx_temp=lwrst_clustMembsCell{k};
            lwrst_dist_temp=lwrst_clust_mat(:,lwrst_clustMembsCell{k});
            lwrst_dist_temp=transpose(lwrst_dist_temp);
            lwrst_dist_mat=pdist(lwrst_dist_temp,'mahalanobis');
            lwrst_dist_mat=squareform(lwrst_dist_mat);
            
            Eps=0.72;
            MinPts=120;

            lwrst_Clust = DBSCAN(lwrst_dist_mat,Eps,MinPts);
            lwrst_indx_clust{k,1}=lwrst_indx_temp(lwrst_Clust~=0);
        end
        
%%        
        Lwrst_indx = horzcat(lwrst_indx_clust{:});
        Rwrst_indx = horzcat(rwrst_indx_clust{:});
        
        wrst_indx=intersect(Lwrst_indx,Rwrst_indx);
                
%% Right Wrist
        
        % Unclean Traj
        figure;
        hold on
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
        for k = 1:min(rwrst_numClust,length(cVec))
            cop_num=strcat("COP",num2str(k));
            cent_num=strcat("Center",num2str(k));
            myMembers = rwrst_clustMembsCell{k};
            myClustCen = rwrst_clustCent(:,k);
            plot(rwrst_clust_mat(1,myMembers),rwrst_clust_mat(2,myMembers),[cVec(k) '.'],'DisplayName',cop_num);
            plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k),'MarkerSize',10,'DisplayName',cent_num);
            legend('-DynamicLegend');
            legend('show');
        end
        hold off;
        title("Unclean RWRIST Trajectory")
        xlabel('Right Wrist Y'),ylabel('Right Wrist Y');

        
        % Unclean Elevation
        figure;
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
        for k = 1:min(rwrst_numClust,length(cVec))
            subplot(2,1,1);
            hold on;
            plot(stand_time(rwrst_clustMembsCell{k}),stand_RWRSTz(rwrst_clustMembsCell{k}),[cVec(k) '.']);
            title("RIGHT WRIST Unclean")
            ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;

            subplot(2,1,2);
            hold on;
            plot(stand_time(rwrst_clustMembsCell{k}),stand_LWRSTz(rwrst_clustMembsCell{k}),[cVec(k) '.']);
            title("LEFT WRIST Unclean")
            ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;
        end
        
        
        % Clean Traj
        figure;
        hold on
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
        for k = 1:min(rwrst_numClust,length(cVec))
            cop_num=strcat("COP",num2str(k));
            cent_num=strcat("Center",num2str(k));
            if length(rwrst_indx_clust{k})>1
                myMembers = rwrst_indx_clust{k};
                temp_clust=rwrst_clust_mat(:,myMembers);
                myClustCen = mean(temp_clust,2);
                plot(rwrst_clust_mat(1,myMembers),rwrst_clust_mat(2,myMembers),[cVec(k) '.'],'DisplayName',cop_num);
                plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k),'MarkerSize',10,'DisplayName',cent_num);
                legend('-DynamicLegend');
                legend('show');
                clear temp_clust
            end
        end
        hold off;
        title("Clean Trajectory")
        xlabel('Right Wrist Y'),ylabel('Right Wrist Y');
        
        % Clean Elevation
        figure;
        for k = 1:min(rwrst_numClust,length(cVec))
            if length(rwrst_indx_clust{k})>1
                subplot(2,1,1);
                hold on;
                plot(stand_time(rwrst_indx_clust{k}),stand_RWRSTz(rwrst_indx_clust{k}),[cVec(k) '.']);
                title("RIGHT WRIST Clean")
                ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;

                subplot(2,1,2);
                hold on;
                plot(stand_time(rwrst_indx_clust{k}),stand_LWRSTz(rwrst_indx_clust{k}),[cVec(k) '.']);
                title("LEFT WRIST Clean")
                ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;
            end
        end
     
        
%% Left Wrist
        
        % Unclean Traj
        figure;
        hold on
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
        for k = 1:min(lwrst_numClust,length(cVec))
            cop_num=strcat("COP",num2str(k));
            cent_num=strcat("Center",num2str(k));
            myMembers = lwrst_clustMembsCell{k};
            myClustCen = lwrst_clustCent(:,k);
            plot(lwrst_clust_mat(1,myMembers),lwrst_clust_mat(2,myMembers),[cVec(k) '.'],'DisplayName',cop_num);
            plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k),'MarkerSize',10,'DisplayName',cent_num);
            legend('-DynamicLegend');
            legend('show');
        end
        hold off;
        title("Unclean LWRIST Trajectory")
        xlabel('Left Wrist Y'),ylabel('Left Wrist Y');
        
        % Unclean Elevation
        figure;
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
        for k = 1:min(lwrst_numClust,length(cVec))
            subplot(2,1,1);
            hold on;
            plot(stand_time(lwrst_clustMembsCell{k}),stand_RWRSTz(lwrst_clustMembsCell{k}),[cVec(k) '.']);
            title("RIGHT WRIST Unclean")
            ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;

            subplot(2,1,2);
            hold on;
            plot(stand_time(lwrst_clustMembsCell{k}),stand_LWRSTz(lwrst_clustMembsCell{k}),[cVec(k) '.']);
            title("LEFT WRIST Unclean")
            ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;
        end
        
        
        % Clean Traj
        figure;
        hold on
        cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
        for k = 1:min(lwrst_numClust,length(cVec))
            cop_num=strcat("COP",num2str(k));
            cent_num=strcat("Center",num2str(k));
            if length(lwrst_indx_clust{k})>1
                myMembers = lwrst_indx_clust{k};
                temp_clust=lwrst_clust_mat(:,myMembers);
                myClustCen = mean(temp_clust,2);
                plot(lwrst_clust_mat(1,myMembers),lwrst_clust_mat(2,myMembers),[cVec(k) '.'],'DisplayName',cop_num);
                plot(myClustCen(1),myClustCen(2),'o','MarkerEdgeColor','k','MarkerFaceColor',cVec(k),'MarkerSize',10,'DisplayName',cent_num);
                legend('-DynamicLegend');
                legend('show');
                clear temp_clust
            end
        end
        hold off;
        title("Clean Trajectory")
        xlabel('Left Wrist Y'),ylabel('Left Wrist Y');
        
        
        % Clean Elevation
        figure;
        for k = 1:min(lwrst_numClust,length(cVec))
            if length(lwrst_indx_clust{k})>1
                subplot(2,1,1);
                hold on;
                plot(stand_time(lwrst_indx_clust{k}),stand_RWRSTz(lwrst_indx_clust{k}),[cVec(k) '.']);
                title("RIGHT WRIST Clean")
                ylabel('Right Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;

                subplot(2,1,2);
                hold on;
                plot(stand_time(lwrst_indx_clust{k}),stand_LWRSTz(lwrst_indx_clust{k}),[cVec(k) '.']);
                title("Left WRIST Clean")
                ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;
            end
        end
        
%% COMBINE CLEAN

        % Clean Elevation
        figure;
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

%% PLOTS
        
        figure;
        subplot(2,1,1);
        hold on;
        plot(stand_time(wrst_indx),stand_RWRSTz(wrst_indx),[cVec(k) '.']);
        title("RIGHT WRIST Clean")
        ylabel('Shoulder XY'),xlabel(strcat('Time (',stand_range,'secs)'));
        hold off;
        subplot(2,1,2);
        hold on;
        plot(stand_time(wrst_indx),stand_LWRSTz(wrst_indx),[cVec(k) '.']);
        title("LEFT WRIST Clean")
        ylabel('Left Wrist Z'),xlabel(strcat('Time (',stand_range,'secs)'));
        hold off;
        
        
        figure;
        subplot(2,1,1);
        hold on;
        plot(stand_time(wrst_indx),stand_SHO_xy_ang(wrst_indx),['b' '.']);
        title("Shoulder XY Angle")
        ylabel('Shoulder XY'),xlabel(strcat('Time (',stand_range,'secs)'));
        hold off;
        subplot(2,1,2);
        hold on;
        plot(stand_time(wrst_indx),stand_SHO_xz_ang(wrst_indx),['b' '.']);
        title("Shoulder XZ Angle")
        ylabel('Shoulder XZ'),xlabel(strcat('Time (',stand_range,'secs)'));
        hold off;
        
        figure;
        subplot(2,1,1);
        plot(stand_time(wrst_indx),stand_head_ang(wrst_indx));
        title("Head Angle")
        ylabel('Head Angle'),xlabel(strcat('Time (',stand_range,'secs)'));
        subplot(2,1,2);
        plot(stand_time(wrst_indx),stand_back_ang(wrst_indx));
        title("Back Angle")
        ylabel('Back Angle'),xlabel(strcat('Time (',stand_range,'secs)'));
        
%         saveas(gcf,strcat(p_num,'_First_step_detect.jpg'));
        

%% SAVING DATA

        standing=struct;
        if numClust>1
        for k=1:numClust
            if length(indx_clust{k})>1
                standing(k).indx=indx_clust{k};
                standing(k).data_pnts=length(indx_clust{k});
                standing(k).time=length(indx_clust{k})/Fs;
                standing(k).clust_cntr_X=clustCent(1,k);
                standing(k).clust_cntr_Y=clustCent(2,k);

                standing(k).COP_X=stand_COP_X(indx_clust{k});
                standing(k).COP_Y=stand_COP_Y(indx_clust{k});
                standing(k).COP_XY=(((standing(k).COP_X).^2)+((standing(k).COP_X).^2)).^0.5;
                standing(k).COP_XY_ratio=(standing(k).COP_X)./(standing(k).COP_Y);

                standing(k).mean_COP_X=mean(standing(k).COP_X);
                standing(k).mean_COP_Y=mean(standing(k).COP_Y);
                standing(k).mean_COP_XY=mean(standing(k).COP_XY);
                standing(k).mean_COP_XY_ratio=mean(standing(k).COP_XY_ratio);

                standing(k).std_COP_X=std(standing(k).COP_X);
                standing(k).std_COP_Y=std(standing(k).COP_Y);
                standing(k).std_COP_XY=std(standing(k).COP_XY);
                standing(k).std_COP_XY_ratio=std(standing(k).COP_XY_ratio);
            end
        end

        stand_parameters=struct;

        stand_parameters.transition=transitions;

        stand_parameters.rel_mean_COP_X=(sum([standing.mean_COP_X].*[standing.data_pnts]))/sum([standing.data_pnts]);
        stand_parameters.rel_mean_COP_Y=(sum([standing.mean_COP_Y].*[standing.data_pnts]))/sum([standing.data_pnts]);

        stand_parameters.dir_mean_COP_X=mean([standing.mean_COP_X]);
        stand_parameters.dir_mean_COP_Y=mean([standing.mean_COP_Y]);

        save(strcat(p_num,'_stand_parameters.mat'),'stand_parameters');
        end

%%
        save(strcat(p_num,'_standing.mat'),'standing');

        fprintf('MAT file save\n');
%         clearvars -except main_workdir exp_workdir subject_num exp_id;
    end
end

      

%%
% COP Vel and root mean square error in anterior-posterior direction
fullCOP_ap_vel=sum((abs(diff(stand_COP_X))./length(stand_COP_X)).*Fs);
fullCOP_ap_rmse=sqrt(sum((stand_COP_X-mean(stand_COP_X)).^2./length(stand_COP_X)));

% COP Vel and root mean square error in medio-lateral direction
fullCOP_ml_vel=sum((abs(diff(stand_COP_Y))./length(stand_COP_Y)).*Fs);
fullCOP_ml_rmse=sqrt(sum((stand_COP_Y-mean(stand_COP_Y)).^2./length(stand_COP_Y)));

