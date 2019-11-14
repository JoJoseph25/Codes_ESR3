clear all; clc; close all;

main_workdir = 'G:\Socrates\Codes'; % Lab Computer  

cd(main_workdir);

exp_workdir='G:\Socrates\Exp_Data\Henrique_data\Exp_data';

% subject_num={'1005';'1006';'1007';'1008';...
%              '1009';'1011';'1012';'1013';...
%              '1014';'1015';'1016';'1017';...
%              '1018';'1019';'1020';'1021';...
%              '1022';'1023';'1024';'1025';...
%              '1026';'1027';'1028'...
%              };



subject_num={'1018'};

exp_id={'dy7'};

exceptions={'1010_dy6','1010_dy5'};

%% PATH CHNG
for sub=1:length(subject_num)
    for exp=1:length(exp_id)
        
        cd(main_workdir);

        x_num=subject_num{sub};

        sub_file=strcat(exp_workdir,'\',x_num);

        addpath(main_workdir)
        addpath(sub_file)

        cd(sub_file)

        p_num=strcat(x_num,'_',exp_id{exp});
        
        p_num
        
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

%% Quiet Standing COP

        % Time
        f30_sec=30*120; % Remove first 30-sec 
        dt = 1/Fs; % Fs=120Hz
        t = FP_data.TIME(:); % Time

        % Time point when the subject moves from Force Plate

        % Value of Force Z around zero (first value Force_Z>-5 + RHEEL\LHEEL>2)
        Force_Z = FP_data.Force_Z(f30_sec:end); % Force Z
        RHEEL_X = Mrk_Data.RHEELX(f30_sec:end); % RHEEL X
        LHEEL_X = Mrk_Data.LHEELX(f30_sec:end); % LHEEL X
        
        indx_frst_stp = find((Force_Z>-5 & (RHEEL_X >10 | LHEEL_X >10)),1,'first'); % Indicator First Step Initiation
        indx_frst_stp = indx_frst_stp + f30_sec;
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

        % COP Data
        COPx = FP_data.COP_X(:); % COP X
        COPx(isnan(COPx))= nanmean(COPx); % Fill the NaN Values ? with the mean of the values of that column that are not NaN
        COPy = FP_data.COP_Y(:); % COP Y
        COPy(isnan(COPy))= nanmean(COPy); % Fill the NaN Values ?

        num_rec = length(COPx); % Number of records 

        % Use the filter
        COPx_Filter = filtfilt(b_fp,a_fp,COPx)/1; % m
        COPy_Filter = filtfilt(b_fp,a_fp,COPy)/1; % m

        % COP XY - for the plot
        COPxy_Filter = (((COPx_Filter.^2)+(COPy_Filter.^2)).^0.5); % m

        % COP Quiet Standing 
        % from 30 seconds - till 5 sec before the subject leaves the FP
        stand_range=num2str((ind_end-ind_start)/Fs);

        stand_time = FP_data.TIME(ind_start:ind_end,:);

        stand_COP_X = COPx_Filter(ind_start:ind_end);
        stand_COP_Y = COPy_Filter(ind_start:ind_end);
        stand_COP_XY= COPxy_Filter(ind_start:ind_end);

        stand_COP_X_vel(1)=0;
        stand_COP_Y_vel(1)=0;
        stand_COP_X_vel(2:length(stand_COP_X))=diff(stand_COP_X)/dt;
        stand_COP_Y_vel(2:length(stand_COP_Y))=diff(stand_COP_Y)/dt;
       
%% PLOTS
        
        figure;
        plotedit('on')
        subplot(3,1,1);
        plot(FP_data.TIME,FP_data.Force_Z);
        line([time_frst_stp time_frst_stp],[-1000 1000],'Color',[1 0 0])
        title("Force Z across Recording")
        ylabel('Force Z-Axis'),xlabel(strcat('Time (',stand_range,'secs)'));
        
        subplot(3,1,2);
        plot(Mrk_Data.Time,Mrk_Data.RHEELX);
        line([time_frst_stp time_frst_stp],[-5000 5000],'Color',[1 0 0])
        title("RHEEL across Recording")
        ylabel('RHEEL X-Axis'),xlabel(strcat('Time (',stand_range,'secs)'));
        
        subplot(3,1,3);
        plot(Mrk_Data.Time,Mrk_Data.LHEELX);
        line([time_frst_stp time_frst_stp],[-5000 5000],'Color',[1 0 0])
        title("LHEEL across Recording")
        ylabel('LHEEL X-Axis'),xlabel(strcat('Time (',stand_range,'secs)'));
        saveas(gcf,strcat(p_num,'_First_step_detect.jpg'));
        close;
        
        figure;
        plot(FP_data.COP_X,FP_data.COP_Y);
        title("Full COP X & Y")
        xlabel('COP X [mm] (medial-lateral)'),ylabel('COP Y [mm] (anterior-posterior)');
        saveas(gcf,strcat(p_num,'_fullCOP.jpg'));
        close;


        figure;
        plot(stand_COP_X,stand_COP_Y);
        title(strcat('COP X & Y 30:-5 sec before First Step (',stand_range,'secs)'));
        xlabel('COP X [mm] (medial-lateral)'),ylabel('COP Y [mm] (anterior-posterior)');
        saveas(gcf,strcat(p_num,'_standCOP.jpg'));
        close;

        figure;
        subplot(3,1,1);
        plot(stand_time,stand_COP_X);
        title("COP X 30:-5 sec before First Step")
        ylabel('COP X [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
        subplot(3,1,2);
        plot(stand_time,stand_COP_Y);
        title("COP X 30:-5 sec before First Step")
        ylabel('COP Y [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
        subplot(3,1,3);
        plot(stand_time,stand_COP_XY);
        title("COP XY 30:-5 sec before First Step")
        ylabel('COP XY [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
        saveas(gcf,strcat(p_num,'_COP_time.jpg'));
        close;

        fprintf('Initial Plot\n');
        
        clear FP_data Mrk_Data LHEEL_X RHEEL_X opts_FP opts_mrk 
        
        clearvars -except p_num main_workdir exp_workdir subject_num exp_id exceptions indx* stand* sub*
%% COP CLustering
        
        
        clust_mat=horzcat(stand_COP_X,stand_COP_Y);

        clust_mat=transpose(clust_mat);

        if ismember(p_num,exceptions)
            min_dist=15;    % Special Cases minimum distance
        else
            min_dist=25;    % Normal minimum distance
        end    

        [clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(clust_mat,min_dist);
                
        numClust = length(clustMembsCell);
        for k=1:numClust
            indx_temp=clustMembsCell{k};
            dist_temp=clust_mat(:,clustMembsCell{k});
            dist_temp=transpose(dist_temp);
            dist_mat=pdist(dist_temp,'mahalanobis');
            clear dist_temp
            dist_mat=squareform(dist_mat);
            
            Eps=0.72;
            MinPts=240;

            Clust = DBSCAN(dist_mat,Eps,MinPts);
            indx_clust{k,1}=indx_temp(Clust~=0);
        end   


        fprintf('Clustering\n');


%% CLUSTER PLOTS

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
        saveas(gcf,strcat(p_num,'_COP_clust.jpg'));
        close;


        figure;
        for k = 1:min(numClust,length(cVec))
            subplot(3,1,1);
            hold on;
            plot(stand_time(clustMembsCell{k}),stand_COP_X(clustMembsCell{k}),[cVec(k) '.']);
            title("COP X 30:-5 sec before First Step")
            ylabel('COP X [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;

            subplot(3,1,2);
            hold on;
            plot(stand_time(clustMembsCell{k}),stand_COP_Y(clustMembsCell{k}),[cVec(k) '.']);
            title("COP Y 30:-5 sec before First Step")
            ylabel('COP Y [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;

            subplot(3,1,3);
            hold on;
            plot(stand_time(clustMembsCell{k}),stand_COP_XY(clustMembsCell{k}),[cVec(k) '.']);
            title("COP XY 30:-5 sec before First Step")
            ylabel('COP XY [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
            hold off;
        end
        saveas(gcf,strcat(p_num,'_COP_clust_time.jpg'));
        close;

        fprintf('Cluster Plots\n');

%% CLEAN CLUSTER PLOTS

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
        saveas(gcf,strcat(p_num,'_clean_COP_clust.jpg'));
        close;


        figure;
        for k = 1:min(numClust,length(cVec))
            if length(indx_clust{k})>1
                subplot(3,1,1);
                hold on;
                plot(stand_time(indx_clust{k}),stand_COP_X(indx_clust{k}),[cVec(k) '.']);
                title("COP X 30:-5 sec before First Step")
                ylabel('COP X [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;

                subplot(3,1,2);
                hold on;
                plot(stand_time(indx_clust{k}),stand_COP_Y(indx_clust{k}),[cVec(k) '.']);
                title("COP Y 30:-5 sec before First Step")
                ylabel('COP Y [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;

                subplot(3,1,3);
                hold on;
                plot(stand_time(indx_clust{k}),stand_COP_XY(indx_clust{k}),[cVec(k) '.']);
                title("COP XY 30:-5 sec before First Step")
                ylabel('COP XY [mm]'),xlabel(strcat('Time (',stand_range,'secs)'));
                hold off;
            end
        end
        saveas(gcf,strcat(p_num,'_clean_COP_clust_time.jpg'));
        close;

        fprintf('Clean Cluster\n');

%% NUM OF TRANSITION

        if numClust>1
            trans_indx=[];
            Clust_id=[];
        
            for i=1:numClust
               
               clear temp_id;
               trans_indx =[trans_indx (indx_clust{i})];
               temp_id= ones(1,length(indx_clust{i}));
               Clust_id =[Clust_id (temp_id.*i)];
    
            end

            clust_trans=[trans_indx;Clust_id]';
            clust_trans_sort=sortrows(clust_trans);

            trans_mat=diff(clust_trans);

            transitions=sum((trans_mat(:,1)>1) | (trans_mat(:,2)>0));
            transitions
        else
            transitions=0;
            
        end

        fprintf('Number of Transitions\n');

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
%}
        fprintf('MAT file save\n');
        clearvars -except main_workdir exp_workdir subject_num exp_id sub exp;
    end
end

      

%%
% COP Vel and root mean square error in anterior-posterior direction
fullCOP_ap_vel=sum((abs(diff(stand_COP_X))./length(stand_COP_X)).*Fs);
fullCOP_ap_rmse=sqrt(sum((stand_COP_X-mean(stand_COP_X)).^2./length(stand_COP_X)));

% COP Vel and root mean square error in medio-lateral direction
fullCOP_ml_vel=sum((abs(diff(stand_COP_Y))./length(stand_COP_Y)).*Fs);
fullCOP_ml_rmse=sqrt(sum((stand_COP_Y-mean(stand_COP_Y)).^2./length(stand_COP_Y)));

