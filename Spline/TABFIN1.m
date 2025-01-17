%AGGIUNGERE IL FOLDER DATI MATLAB E TUTTI I SUBFOLDERS AL PATH
addpath(genpath('D:\SIMONEUTENTE\DESKTOP\PROG_QUA\ESPERIMENTO\QUASISTATIC\RISULTATI_MATLAB'))
clear
close all
clc
directory = dir ('D:\SIMONEUTENTE\DESKTOP\PROG_QUA\ESPERIMENTO\QUASISTATIC\RISULTATI_MATLAB'); % bisogna stare in questa directory
iii=size(directory,1);


TABFIN=cell(3250,569);
TABFIN(1,:)={'ID','SubjectHeight','SubjectWeight','Sex','Vo2Max','BoxMass','LiftingHeight','Loweringheight',...
    'LIF_Peak_Back_MomX','LIF_Peak_Lsho_MomX','LIF_Peak_Rsho_MomX','LIF_Peak_Lelb_MomX','LIF_Peak_Relb_MomX','LIF_Peak_Lwrist_MomX','LIF_Peak_Rwrist_MomX',...
    'LIF_Peak_Back_MomNET','LIF_Peak_Lsho_MomNET','LIF_Peak_Rsho_MomNET','LIF_Peak_Lelb_MomNET','LIF_Peak_Relb_MomNET','LIF_Peak_Lwrist_MomNET','LIF_Peak_Rwrist_MomNET',...
    'LIF_Cumu_Back_MomX','LIF_Cumu_Lsho_MomX','LIF_Cumu_Rsho_MomX','LIF_Cumu_Lelb_MomX','LIF_Cumu_Relb_MomX','LIF_Cumu_Lwrist_MomX','LIF_Cumu_Rwrist_MomX',...
    'LIF_Cumu_Back_MomNET','LIF_Cumu_Lsho_MomNET','LIF_Cumu_Rsho_MomNET','LIF_Cumu_Lelb_MomNET','LIF_Cumu_Relb_MomNET','LIF_Cumu_Lwrist_MomNET','LIF_Cumu_Rwrist_MomNET',...
    'LIF_Peak_Trunk_Angle','LIF_Peak_Lsho_Angle','LIF_Peak_Rsho_Angle','LIF_Peak_Lelb_Angle','LIF_Peak_Relb_Angle','LIF_Peak_Lwrist_Angle','LIF_Peak_Rwrist_Angle',...
    'LIF_Peak_Lhip_Angle','LIF_Peak_Rhip_Angle','LIF_Peak_Lknee_Angle','LIF_Peak_Rknee_Angle','LIF_Peak_Lankle_Angle','LIF_Peak_Rankle_Angle',...
    'LIF_Peak_Trunk_Vel','LIF_Peak_Lsho_Vel','LIF_Peak_Rsho_Vel','LIF_Peak_Lelb_Vel','LIF_Peak_Relb_Vel','LIF_Peak_Lwrist_Vel','LIF_Peak_Rwrist_Vel',...
    'LIF_Peak_Lhip_Vel','LIF_Peak_Rhip_Vel','LIF_Peak_Lknee_Vel','LIF_Peak_Rknee_Vel','LIF_Peak_Lankle_Vel','LIF_Peak_Rankle_Vel',...
    'LIF_Peak_Trunk_Acc','LIF_Peak_Lsho_Acc','LIF_Peak_Rsho_Acc','LIF_Peak_Lelb_Acc','LIF_Peak_Relb_Acc','LIF_Peak_Lwrist_Acc','LIF_Peak_Rwrist_Acc',...
    'LIF_Peak_Lhip_Acc','LIF_Peak_Rhip_Acc','LIF_Peak_Lknee_Acc','LIF_Peak_Rknee_Acc','LIF_Peak_Lankle_Acc','LIF_Peak_Rankle_Acc',...
    'LIF_PMB_Trunk_Angle','LIF_PMB_Lsho_Angle','LIF_PMB_Rsho_Angle','LIF_PMB_Lelb_Angle','LIF_PMB_Relb_Angle','LIF_PMB_Lwrist_Angle','LIF_PMB_Rwrist_Angle',...
    'LIF_PMB_Lhip_Angle','LIF_PMB_Rhip_Angle','LIF_PMB_Lknee_Angle','LIF_PMB_Rknee_Angle','LIF_PMB_Lankle_Angle','LIF_PMB_Rankle_Angle',...
    'LIF_PMLSH_Trunk_Angle','LIF_PMLSH_Lsho_Angle','LIF_PMLSH_Rsho_Angle','LIF_PMLSH_Lelb_Angle','LIF_PMLSH_Relb_Angle','LIF_PMLSH_Lwrist_Angle','LIF_PMLSH_Rwrist_Angle',...
    'LIF_PMLSH_Lhip_Angle','LIF_PMLSH_Rhip_Angle','LIF_PMLSH_Lknee_Angle','LIF_PMLSH_Rknee_Angle','LIF_PMLSH_Lankle_Angle','LIF_PMLSH_Rankle_Angle',...
    'LIF_PMRSH_Trunk_Angle','LIF_PMRSH_Lsho_Angle','LIF_PMRSH_Rsho_Angle','LIF_PMRSH_Lelb_Angle','LIF_PMRSH_Relb_Angle','LIF_PMRSH_Lwrist_Angle','LIF_PMRSH_Rwrist_Angle',...
    'LIF_PMRSH_Lhip_Angle','LIF_PMRSH_Rhip_Angle','LIF_PMRSH_Lknee_Angle','LIF_PMRSH_Rknee_Angle','LIF_PMRSH_Lankle_Angle','LIF_PMRSH_Rankle_Angle',...
    'LIF_PMB_Trunk_Vel','LIF_PMB_Lsho_Vel','LIF_PMB_Rsho_Vel','LIF_PMB_Lelb_Vel','LIF_PMB_Relb_Vel','LIF_PMB_Lwrist_Vel','LIF_PMB_Rwrist_Vel',...
    'LIF_PMB_Lhip_Vel','LIF_PMB_Rhip_Vel','LIF_PMB_Lknee_Vel','LIF_PMB_Rknee_Vel','LIF_PMB_Lankle_Vel','LIF_PMB_Rankle_Vel',...
    'LIF_PMLSH_Trunk_Vel','LIF_PMLSH_Lsho_Vel','LIF_PMLSH_Rsho_Vel','LIF_PMLSH_Lelb_Vel','LIF_PMLSH_Relb_Vel','LIF_PMLSH_Lwrist_Vel','LIF_PMLSH_Rwrist_Vel',...
    'LIF_PMLSH_Lhip_Vel','LIF_PMLSH_Rhip_Vel','LIF_PMLSH_Lknee_Vel','LIF_PMLSH_Rknee_Vel','LIF_PMLSH_Lankle_Vel','LIF_PMLSH_Rankle_Vel',...
    'LIF_PMRSH_Trunk_Vel','LIF_PMRSH_Lsho_Vel','LIF_PMRSH_Rsho_Vel','LIF_PMRSH_Lelb_Vel','LIF_PMRSH_Relb_Vel','LIF_PMRSH_Lwrist_Vel','LIF_PMRSH_Rwrist_Vel',...
    'LIF_PMRSH_Lhip_Vel','LIF_PMRSH_Rhip_Vel','LIF_PMRSH_Lknee_Vel','LIF_PMRSH_Rknee_Vel','LIF_PMRSH_Lankle_Vel','LIF_PMRSH_Rankle_Vel',...
    'LIF_PMB_Trunk_Acc','LIF_PMB_Lsho_Acc','LIF_PMB_Rsho_Acc','LIF_PMB_Lelb_Acc','LIF_PMB_Relb_Acc','LIF_PMB_Lwrist_Acc','LIF_PMB_Rwrist_Acc',...
    'LIF_PMB_Lhip_Acc','LIF_PMB_Rhip_Acc','LIF_PMB_Lknee_Acc','LIF_PMB_Rknee_Acc','LIF_PMB_Lankle_Acc','LIF_PMB_Rankle_Acc',...
    'LIF_PMLSH_Trunk_Acc','LIF_PMLSH_Lsho_Acc','LIF_PMLSH_Rsho_Acc','LIF_PMLSH_Lelb_Acc','LIF_PMLSH_Relb_Acc','LIF_PMLSH_Lwrist_Acc','LIF_PMLSH_Rwrist_Acc',...
    'LIF_PMLSH_Lhip_Acc','LIF_PMLSH_Rhip_Acc','LIF_PMLSH_Lknee_Acc','LIF_PMLSH_Rknee_Acc','LIF_PMLSH_Lankle_Acc','LIF_PMLSH_Rankle_Acc',...
    'LIF_PMRSH_Trunk_Acc','LIF_PMRSH_Lsho_Acc','LIF_PMRSH_Rsho_Acc','LIF_PMRSH_Lelb_Acc','LIF_PMRSH_Relb_Acc','LIF_PMRSH_Lwrist_Acc','LIF_PMRSH_Rwrist_Acc',...
    'LIF_PMRSH_Lhip_Acc','LIF_PMRSH_Rhip_Acc','LIF_PMRSH_Lknee_Acc','LIF_PMRSH_Rknee_Acc','LIF_PMRSH_Lankle_Acc','LIF_PMRSH_Rankle_Acc',...
    'LIF_Peak_Arm','LIF_PM_Arm',...
    'CAR_Peak_Back_MomX','CAR_Peak_Lsho_MomX','CAR_Peak_Rsho_MomX','CAR_Peak_Lelb_MomX','CAR_Peak_Relb_MomX','CAR_Peak_Lwrist_MomX','CAR_Peak_Rwrist_MomX',...
    'CAR_Peak_Back_MomNET','CAR_Peak_Lsho_MomNET','CAR_Peak_Rsho_MomNET','CAR_Peak_Lelb_MomNET','CAR_Peak_Relb_MomNET','CAR_Peak_Lwrist_MomNET','CAR_Peak_Rwrist_MomNET',...
    'CAR_Cumu_Back_MomX','CAR_Cumu_Lsho_MomX','CAR_Cumu_Rsho_MomX','CAR_Cumu_Lelb_MomX','CAR_Cumu_Relb_MomX','CAR_Cumu_Lwrist_MomX','CAR_Cumu_Rwrist_MomX',...
    'CAR_Cumu_Back_MomNET','CAR_Cumu_Lsho_MomNET','CAR_Cumu_Rsho_MomNET','CAR_Cumu_Lelb_MomNET','CAR_Cumu_Relb_MomNET','CAR_Cumu_Lwrist_MomNET','CAR_Cumu_Rwrist_MomNET',...
    'CAR_Peak_Trunk_Angle','CAR_Peak_Lsho_Angle','CAR_Peak_Rsho_Angle','CAR_Peak_Lelb_Angle','CAR_Peak_Relb_Angle','CAR_Peak_Lwrist_Angle','CAR_Peak_Rwrist_Angle',...
    'CAR_Peak_Lhip_Angle','CAR_Peak_Rhip_Angle','CAR_Peak_Lknee_Angle','CAR_Peak_Rknee_Angle','CAR_Peak_Lankle_Angle','CAR_Peak_Rankle_Angle',...
    'CAR_Peak_Trunk_Vel','CAR_Peak_Lsho_Vel','CAR_Peak_Rsho_Vel','CAR_Peak_Lelb_Vel','CAR_Peak_Relb_Vel','CAR_Peak_Lwrist_Vel','CAR_Peak_Rwrist_Vel',...
    'CAR_Peak_Lhip_Vel','CAR_Peak_Rhip_Vel','CAR_Peak_Lknee_Vel','CAR_Peak_Rknee_Vel','CAR_Peak_Lankle_Vel','CAR_Peak_Rankle_Vel',...
    'CAR_Peak_Trunk_Acc','CAR_Peak_Lsho_Acc','CAR_Peak_Rsho_Acc','CAR_Peak_Lelb_Acc','CAR_Peak_Relb_Acc','CAR_Peak_Lwrist_Acc','CAR_Peak_Rwrist_Acc',...
    'CAR_Peak_Lhip_Acc','CAR_Peak_Rhip_Acc','CAR_Peak_Lknee_Acc','CAR_Peak_Rknee_Acc','CAR_Peak_Lankle_Acc','CAR_Peak_Rankle_Acc',...
    'CAR_PMB_Trunk_Angle','CAR_PMB_Lsho_Angle','CAR_PMB_Rsho_Angle','CAR_PMB_Lelb_Angle','CAR_PMB_Relb_Angle','CAR_PMB_Lwrist_Angle','CAR_PMB_Rwrist_Angle',...
    'CAR_PMB_Lhip_Angle','CAR_PMB_Rhip_Angle','CAR_PMB_Lknee_Angle','CAR_PMB_Rknee_Angle','CAR_PMB_Lankle_Angle','CAR_PMB_Rankle_Angle',...
    'CAR_PMLSH_Trunk_Angle','CAR_PMLSH_Lsho_Angle','CAR_PMLSH_Rsho_Angle','CAR_PMLSH_Lelb_Angle','CAR_PMLSH_Relb_Angle','CAR_PMLSH_Lwrist_Angle','CAR_PMLSH_Rwrist_Angle',...
    'CAR_PMLSH_Lhip_Angle','CAR_PMLSH_Rhip_Angle','CAR_PMLSH_Lknee_Angle','CAR_PMLSH_Rknee_Angle','CAR_PMLSH_Lankle_Angle','CAR_PMLSH_Rankle_Angle',...
    'CAR_PMRSH_Trunk_Angle','CAR_PMRSH_Lsho_Angle','CAR_PMRSH_Rsho_Angle','CAR_PMRSH_Lelb_Angle','CAR_PMRSH_Relb_Angle','CAR_PMRSH_Lwrist_Angle','CAR_PMRSH_Rwrist_Angle',...
    'CAR_PMRSH_Lhip_Angle','CAR_PMRSH_Rhip_Angle','CAR_PMRSH_Lknee_Angle','CAR_PMRSH_Rknee_Angle','CAR_PMRSH_Lankle_Angle','CAR_PMRSH_Rankle_Angle',...
    'CAR_PMB_Trunk_Vel','CAR_PMB_Lsho_Vel','CAR_PMB_Rsho_Vel','CAR_PMB_Lelb_Vel','CAR_PMB_Relb_Vel','CAR_PMB_Lwrist_Vel','CAR_PMB_Rwrist_Vel',...
    'CAR_PMB_Lhip_Vel','CAR_PMB_Rhip_Vel','CAR_PMB_Lknee_Vel','CAR_PMB_Rknee_Vel','CAR_PMB_Lankle_Vel','CAR_PMB_Rankle_Vel',...
    'CAR_PMLSH_Trunk_Vel','CAR_PMLSH_Lsho_Vel','CAR_PMLSH_Rsho_Vel','CAR_PMLSH_Lelb_Vel','CAR_PMLSH_Relb_Vel','CAR_PMLSH_Lwrist_Vel','CAR_PMLSH_Rwrist_Vel',...
    'CAR_PMLSH_Lhip_Vel','CAR_PMLSH_Rhip_Vel','CAR_PMLSH_Lknee_Vel','CAR_PMLSH_Rknee_Vel','CAR_PMLSH_Lankle_Vel','CAR_PMLSH_Rankle_Vel',...
    'CAR_PMRSH_Trunk_Vel','CAR_PMRSH_Lsho_Vel','CAR_PMRSH_Rsho_Vel','CAR_PMRSH_Lelb_Vel','CAR_PMRSH_Relb_Vel','CAR_PMRSH_Lwrist_Vel','CAR_PMRSH_Rwrist_Vel',...
    'CAR_PMRSH_Lhip_Vel','CAR_PMRSH_Rhip_Vel','CAR_PMRSH_Lknee_Vel','CAR_PMRSH_Rknee_Vel','CAR_PMRSH_Lankle_Vel','CAR_PMRSH_Rankle_Vel',...
    'CAR_PMB_Trunk_Acc','CAR_PMB_Lsho_Acc','CAR_PMB_Rsho_Acc','CAR_PMB_Lelb_Acc','CAR_PMB_Relb_Acc','CAR_PMB_Lwrist_Acc','CAR_PMB_Rwrist_Acc',...
    'CAR_PMB_Lhip_Acc','CAR_PMB_Rhip_Acc','CAR_PMB_Lknee_Acc','CAR_PMB_Rknee_Acc','CAR_PMB_Lankle_Acc','CAR_PMB_Rankle_Acc',...
    'CAR_PMLSH_Trunk_Acc','CAR_PMLSH_Lsho_Acc','CAR_PMLSH_Rsho_Acc','CAR_PMLSH_Lelb_Acc','CAR_PMLSH_Relb_Acc','CAR_PMLSH_Lwrist_Acc','CAR_PMLSH_Rwrist_Acc',...
    'CAR_PMLSH_Lhip_Acc','CAR_PMLSH_Rhip_Acc','CAR_PMLSH_Lknee_Acc','CAR_PMLSH_Rknee_Acc','CAR_PMLSH_Lankle_Acc','CAR_PMLSH_Rankle_Acc',...
    'CAR_PMRSH_Trunk_Acc','CAR_PMRSH_Lsho_Acc','CAR_PMRSH_Rsho_Acc','CAR_PMRSH_Lelb_Acc','CAR_PMRSH_Relb_Acc','CAR_PMRSH_Lwrist_Acc','CAR_PMRSH_Rwrist_Acc',...
    'CAR_PMRSH_Lhip_Acc','CAR_PMRSH_Rhip_Acc','CAR_PMRSH_Lknee_Acc','CAR_PMRSH_Rknee_Acc','CAR_PMRSH_Lankle_Acc','CAR_PMRSH_Rankle_Acc',...
    'CAR_Peak_Arm','CAR_PM_Arm',...
    'LOW_Peak_Back_MomX','LOW_Peak_Lsho_MomX','LOW_Peak_Rsho_MomX','LOW_Peak_Lelb_MomX','LOW_Peak_Relb_MomX','LOW_Peak_Lwrist_MomX','LOW_Peak_Rwrist_MomX',...
    'LOW_Peak_Back_MomNET','LOW_Peak_Lsho_MomNET','LOW_Peak_Rsho_MomNET','LOW_Peak_Lelb_MomNET','LOW_Peak_Relb_MomNET','LOW_Peak_Lwrist_MomNET','LOW_Peak_Rwrist_MomNET',...
    'LOW_Cumu_Back_MomX','LOW_Cumu_Lsho_MomX','LOW_Cumu_Rsho_MomX','LOW_Cumu_Lelb_MomX','LOW_Cumu_Relb_MomX','LOW_Cumu_Lwrist_MomX','LOW_Cumu_Rwrist_MomX',...
    'LOW_Cumu_Back_MomNET','LOW_Cumu_Lsho_MomNET','LOW_Cumu_Rsho_MomNET','LOW_Cumu_Lelb_MomNET','LOW_Cumu_Relb_MomNET','LOW_Cumu_Lwrist_MomNET','LOW_Cumu_Rwrist_MomNET',...
    'LOW_Peak_Trunk_Angle','LOW_Peak_Lsho_Angle','LOW_Peak_Rsho_Angle','LOW_Peak_Lelb_Angle','LOW_Peak_Relb_Angle','LOW_Peak_Lwrist_Angle','LOW_Peak_Rwrist_Angle',...
    'LOW_Peak_Lhip_Angle','LOW_Peak_Rhip_Angle','LOW_Peak_Lknee_Angle','LOW_Peak_Rknee_Angle','LOW_Peak_Lankle_Angle','LOW_Peak_Rankle_Angle',...
    'LOW_Peak_Trunk_Vel','LOW_Peak_Lsho_Vel','LOW_Peak_Rsho_Vel','LOW_Peak_Lelb_Vel','LOW_Peak_Relb_Vel','LOW_Peak_Lwrist_Vel','LOW_Peak_Rwrist_Vel',...
    'LOW_Peak_Lhip_Vel','LOW_Peak_Rhip_Vel','LOW_Peak_Lknee_Vel','LOW_Peak_Rknee_Vel','LOW_Peak_Lankle_Vel','LOW_Peak_Rankle_Vel',...
    'LOW_Peak_Trunk_Acc','LOW_Peak_Lsho_Acc','LOW_Peak_Rsho_Acc','LOW_Peak_Lelb_Acc','LOW_Peak_Relb_Acc','LOW_Peak_Lwrist_Acc','LOW_Peak_Rwrist_Acc',...
    'LOW_Peak_Lhip_Acc','LOW_Peak_Rhip_Acc','LOW_Peak_Lknee_Acc','LOW_Peak_Rknee_Acc','LOW_Peak_Lankle_Acc','LOW_Peak_Rankle_Acc',...
    'LOW_PMB_Trunk_Angle','LOW_PMB_Lsho_Angle','LOW_PMB_Rsho_Angle','LOW_PMB_Lelb_Angle','LOW_PMB_Relb_Angle','LOW_PMB_Lwrist_Angle','LOW_PMB_Rwrist_Angle',...
    'LOW_PMB_Lhip_Angle','LOW_PMB_Rhip_Angle','LOW_PMB_Lknee_Angle','LOW_PMB_Rknee_Angle','LOW_PMB_Lankle_Angle','LOW_PMB_Rankle_Angle',...
    'LOW_PMLSH_Trunk_Angle','LOW_PMLSH_Lsho_Angle','LOW_PMLSH_Rsho_Angle','LOW_PMLSH_Lelb_Angle','LOW_PMLSH_Relb_Angle','LOW_PMLSH_Lwrist_Angle','LOW_PMLSH_Rwrist_Angle',...
    'LOW_PMLSH_Lhip_Angle','LOW_PMLSH_Rhip_Angle','LOW_PMLSH_Lknee_Angle','LOW_PMLSH_Rknee_Angle','LOW_PMLSH_Lankle_Angle','LOW_PMLSH_Rankle_Angle',...
    'LOW_PMRSH_Trunk_Angle','LOW_PMRSH_Lsho_Angle','LOW_PMRSH_Rsho_Angle','LOW_PMRSH_Lelb_Angle','LOW_PMRSH_Relb_Angle','LOW_PMRSH_Lwrist_Angle','LOW_PMRSH_Rwrist_Angle',...
    'LOW_PMRSH_Lhip_Angle','LOW_PMRSH_Rhip_Angle','LOW_PMRSH_Lknee_Angle','LOW_PMRSH_Rknee_Angle','LOW_PMRSH_Lankle_Angle','LOW_PMRSH_Rankle_Angle',...
    'LOW_PMB_Trunk_Vel','LOW_PMB_Lsho_Vel','LOW_PMB_Rsho_Vel','LOW_PMB_Lelb_Vel','LOW_PMB_Relb_Vel','LOW_PMB_Lwrist_Vel','LOW_PMB_Rwrist_Vel',...
    'LOW_PMB_Lhip_Vel','LOW_PMB_Rhip_Vel','LOW_PMB_Lknee_Vel','LOW_PMB_Rknee_Vel','LOW_PMB_Lankle_Vel','LOW_PMB_Rankle_Vel',...
    'LOW_PMLSH_Trunk_Vel','LOW_PMLSH_Lsho_Vel','LOW_PMLSH_Rsho_Vel','LOW_PMLSH_Lelb_Vel','LOW_PMLSH_Relb_Vel','LOW_PMLSH_Lwrist_Vel','LOW_PMLSH_Rwrist_Vel',...
    'LOW_PMLSH_Lhip_Vel','LOW_PMLSH_Rhip_Vel','LOW_PMLSH_Lknee_Vel','LOW_PMLSH_Rknee_Vel','LOW_PMLSH_Lankle_Vel','LOW_PMLSH_Rankle_Vel',...
    'LOW_PMRSH_Trunk_Vel','LOW_PMRSH_Lsho_Vel','LOW_PMRSH_Rsho_Vel','LOW_PMRSH_Lelb_Vel','LOW_PMRSH_Relb_Vel','LOW_PMRSH_Lwrist_Vel','LOW_PMRSH_Rwrist_Vel',...
    'LOW_PMRSH_Lhip_Vel','LOW_PMRSH_Rhip_Vel','LOW_PMRSH_Lknee_Vel','LOW_PMRSH_Rknee_Vel','LOW_PMRSH_Lankle_Vel','LOW_PMRSH_Rankle_Vel',...
    'LOW_PMB_Trunk_Acc','LOW_PMB_Lsho_Acc','LOW_PMB_Rsho_Acc','LOW_PMB_Lelb_Acc','LOW_PMB_Relb_Acc','LOW_PMB_Lwrist_Acc','LOW_PMB_Rwrist_Acc',...
    'LOW_PMB_Lhip_Acc','LOW_PMB_Rhip_Acc','LOW_PMB_Lknee_Acc','LOW_PMB_Rknee_Acc','LOW_PMB_Lankle_Acc','LOW_PMB_Rankle_Acc',...
    'LOW_PMLSH_Trunk_Acc','LOW_PMLSH_Lsho_Acc','LOW_PMLSH_Rsho_Acc','LOW_PMLSH_Lelb_Acc','LOW_PMLSH_Relb_Acc','LOW_PMLSH_Lwrist_Acc','LOW_PMLSH_Rwrist_Acc',...
    'LOW_PMLSH_Lhip_Acc','LOW_PMLSH_Rhip_Acc','LOW_PMLSH_Lknee_Acc','LOW_PMLSH_Rknee_Acc','LOW_PMLSH_Lankle_Acc','LOW_PMLSH_Rankle_Acc',...
    'LOW_PMRSH_Trunk_Acc','LOW_PMRSH_Lsho_Acc','LOW_PMRSH_Rsho_Acc','LOW_PMRSH_Lelb_Acc','LOW_PMRSH_Relb_Acc','LOW_PMRSH_Lwrist_Acc','LOW_PMRSH_Rwrist_Acc',...
    'LOW_PMRSH_Lhip_Acc','LOW_PMRSH_Rhip_Acc','LOW_PMRSH_Lknee_Acc','LOW_PMRSH_Rknee_Acc','LOW_PMRSH_Lankle_Acc','LOW_PMRSH_Rankle_Acc',...
    'LOW_Peak_Arm','LOW_PM_Arm',...
    'Lifting_Time','Carrying_time','Lowering_time'};




    instat=2; %indice di riempimento per gli static
    inmot=2; %indice per i motion

for i=3:iii %dentro il FOLDER risultatiMATLAB, ogni elemento � la cart X...X
    sub = char( string(directory(i).folder)+"\"+string(directory(i).name) );
    subdirec = dir (sub);
    qqq=size(subdirec,1);
    
    
    
    
    
    for qq=3:qqq  %dentro il FOLDER X...X, ogni elemento � un tipo Aviv8RIS che contiene un 'datatot'
        load(subdirec(qq).name); %PER OGNI datatot NELLA CARTELLA --> static
        
        
        %ESTRARRE LE INFORMAZIONI COMUNI AD UNO STATIC
        nmot=length(cell2mat(datatot(:,3)));
        for j=1:nmot
            TABFIN{instat,1}= datatot{3,5};   %ID
            TABFIN{instat,2}= datatot{7,5};   %sub height
            TABFIN{instat,3}= datatot{9,5};   %sub weight
            TABFIN{instat,4}= datatot{5,5};   %sex
            TABFIN{instat,6}= datatot{11,5};  %boxmass
            
            TABFIN{instat,7}= datatot{j,3};   %scaff in
            TABFIN{instat,8}= datatot{j,4};  %scaff fin
            instat=instat+1;
        end
        
        
        
        
        %%
        for u=1:nmot% PER OGNI MOTION NELLO STATIC
            j=9;
            
            for h=2:5
                for c=2:8
                    TABFIN{inmot,j}=datatot{u,1}{h,c};
                    j=j+1;
                end
            end
            
            for h=6:17
                for c=2:14
                    TABFIN{inmot,j}=datatot{u,1}{h,c};
                    j=j+1;
                end
            end
            
            TABFIN{inmot,j}=datatot{u,1}{18,15};
            TABFIN{inmot,j+1}=datatot{u,1}{19,15};
            j=j+2;
            
            
            for h=20:23
                for c=2:8
                    TABFIN{inmot,j}=datatot{u,1}{h,c};
                    j=j+1;
                end
            end
            
            
            for h=24:35
                for c=2:14
                    TABFIN{inmot,j}=datatot{u,1}{h,c};
                    j=j+1;
                end
            end
            
            TABFIN{inmot,j}=datatot{u,1}{36,15};
            TABFIN{inmot,j+1}=datatot{u,1}{37,15};
            j=j+2;
            
            
            for h=38:41
                for c=2:8
                    TABFIN{inmot,j}=datatot{u,1}{h,c};
                    j=j+1;
                end
            end
            
            
            for h=42:53
                for c=2:14
                    TABFIN{inmot,j}=datatot{u,1}{h,c};
                    j=j+1;
                end
            end
            
            
            TABFIN{inmot,j}=datatot{u,1}{54,15};
            TABFIN{inmot,j+1}=datatot{u,1}{55,15};
            j=j+2;
            
            
            
            TABFIN{inmot,j}=datatot{u,6};
            TABFIN{inmot,j+1}=datatot{u,7};
            TABFIN{inmot,j+2}=datatot{u,8};
            
            
            inmot=inmot+1;
            
            
        end   %U end-for dei motion files
        
        
       
        
        
        clearvars -except directory iii TABFIN i sub subdirec qqq qq instat inmot
        
    end
    
end



path=string(subdirec(qq).folder);
pars=strsplit(path,'RISULTATI_MATLAB');
path =pars(1,1);
TABLEQUASIS=TABFIN;
path1= path + "TABLEQUASIS";
path1=char(path1);
save(path1, 'TABLEQUASIS')


path2= path + "TABLEQUASIS.xlsx";
path2 = char(path2) ;
xlswrite(path2,TABLEQUASIS)


