%AGGIUNGERE IL FOLDER DATI MATLAB E TUTTI I SUBFOLDERS AL PATH
addpath(genpath('D:\SIMONEUTENTE\DESKTOP\PROG_QUA\ESPERIMENTO\QUASISTATIC\DATI_MATLAB_INTERPOLATI'))
clear
close all
clc
directory = dir ('D:\SIMONEUTENTE\DESKTOP\PROG_QUA\ESPERIMENTO\QUASISTATIC\DATI_MATLAB_INTERPOLATI'); % bisogna stare in questa directory 
iii=size(directory,1);

%mettere in colonna i nomi secondo l'ordine dell'id
nametable={'NAME','ID','SEX';'Amit',7,1;  'Aviv',20,1;  'Ayelet',13,0;  'Chen',17,0;  'Eitan',14,1; 'Gil',12,1;   'Gilad',16,1;  'Ido',2,1;   'Inbal',18,0;...
          'Moria',15,0;  'Noy',11,0;  'Orben',19,0;  'Orkaz',3,0;  'Oz',10,1;   'Ravit',1,0; 'Roy',4,1;    'Saar',9,1;  'Shahar',8,0};
%%%%%%%%%%%%%%%%%%%%%%%%%%%
scafftable={'START_BOX_HEIGHT','END_BOX_HEIGHT';20,20;50,50;80,80;110,110;140,140;170,170};
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=3:iii %dentro il FOLDER DATIMATLAB
sub = char( string(directory(i).folder)+"\"+string(directory(i).name) );
subdirec = dir (sub);  
qqq=size(subdirec,1);
    

%%%%%%%%%%%identificare nome e ID%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name_sub=string(directory(i).name);
parsub=strsplit(name_sub,'X');
name_sub =parsub(1,2);
name_sub=char(name_sub);
id_sub = nametable{strcmp(nametable, name_sub),2};
sex_sub=nametable{strcmp(nametable, name_sub),3};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 for qq=3:qqq  %dentro il FOLDER X...X
 load(subdirec(qq).name);%PER OGNI STATIC NELLA CARTELLA
 le0=length(Back_acc);
 motion = cell(le0,1);
 totframes=zeros(le0,1);
 for j=1:le0
  motion{j,1}= ['motion' num2str(j)];
  totframes(j,1)= size(Back_acc{j,1},1);
 end
 datatot=cell(le0,1);
 
 %%%%%%%%%%%%%%%%%%%%%%%CALCOLO IL VALORE REALE DELLA MASSA SCATOLA
 BOXMASS=BOXMASS{1,1}*2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 
 %%
 for u=1:le0% PER OGNI MOTION NELLO STATIC
     
  %%%%VERIFICO I BACKUP%%%%%%%   
  %lift=[START_LIFTING{u},END_LIFTING{u}];carry=[END_LIFTING{u}, START_LOWERING{u}];lower=[START_LOWERING{u}, END_LOWERING{u}];
  time= {START_LIFTING{u}, END_LIFTING{u}, START_LOWERING{u}, END_LOWERING{u}};
  backup= [CROSS_LAB_0{u},CROSS_LAB_1{u}-60, CROSS_LAB_1{u}+60, CROSS_LAB_2{u}];  % backup
  for j=1:length(time)
      if  isempty(time{j})
          time{j}=backup(j);
      end
  end    
 
  if ( time{1}==CROSS_LAB_0{u} || time{2}==time{3} || time{4}==CROSS_LAB_2{u} )
     datatot{u,2}={'Backup used'};
  end
  time=cell2mat(time);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%COSTRUISCO TABELLA%%%%%%%%%%%%%%%%%%%%%%
  datatot{u,1}= cell(55,15);
  datatot{u,1}(1,2:end)={'Back','L_sh','R_sh','L_el','R_el','L_wr','R_wr','L_hip','R_hip','L_kn','R_kn','L_ank','R_ank','ARM' };
  
  datatot{u,1}(2:end,1)={'peak_X_mom_LIFT';'peak_TOT_mom_LIFT';'cum_X_mom_LIFT';'cum_TOT_mom_LIFT';'peak_X_ang_LIFT';'peak_X_vel_LIFT';'peak_X_acc_LIFT';...
    'PM_Back_angle_LIFT'; 'PM_Lsh_angle_LIFT'; 'PM_Rsh_angle_LIFT';   'PM_Back_vel_LIFT'; 'PM_Lsh_vel_LIFT'; 'PM_Rsh_vel_LIFT';  'PM_Back_acc_LIFT'; 'PM_Lsh_acc_LIFT'; 'PM_Rsh_acc_LIFT';
    'Peak_ARM_LIFT';      'PM_Back_ARM_LIFT';
    
    'peak_X_mom_CARRY';'peak_TOT_mom_CARRY';'cum_X_mom_CARRY';'cum_TOT_mom_CARRY';'peak_X_ang_CARRY'; 'peak_X_vel_CARRY';'peak_X_acc_CARRY';
    'PM_Back_angle_CARRY'; 'PM_Lsh_angle_CARRY'; 'PM_Rsh_angle_CARRY';  'PM_Back_vel_CARRY'; 'PM_Lsh_vel_CARRY'; 'PM_Rsh_vel_CARRY';'PM_Back_acc_CARRY'; 'PM_Lsh_acc_CARRY'; 'PM_Rsh_acc_CARRY';
    'Peak_ARM_CARRY';      'PM_Back_ARM_CARRY'; 
    
    'peak_X_mom_LOW';'peak_TOT_mom_LOW';'cum_X_mom_LOW';'cum_TOT_mom_LOW';'peak_X_ang_LOW';'peak_X_vel_LOW';'peak_X_acc_LOW';
    'PM_Back_angle_LOW'; 'PM_Lsh_angle_LOW'; 'PM_Rsh_angle_LOW';  'PM_Back_vel_LOW'; 'PM_Lsh_vel_LOW'; 'PM_Rsh_vel_LOW';'PM_Back_acc_LOW'; 'PM_Lsh_acc_LOW'; 'PM_Rsh_acc_LOW';
    'Peak_ARM_LOW';      'PM_Back_ARM_LOW'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
 
  
  
  for h=1:length(time)-1  % FASI CARRY
  
  %%%%%%%%%%%%%%%%%MOMENTO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  mom= [ abs( Back_moment{u,1}(time(h):time(h+1),:) ),   abs( L_shoulder_moment{u,1}(time(h):time(h+1),:) ),  abs(R_shoulder_moment{u,1}(time(h):time(h+1),:) ),...
       abs(L_elbow_moment{u,1}(time(h):time(h+1),:) ), abs(R_elbow_moment{u,1}(time(h):time(h+1),:) ), abs(L_wrist_moment{u,1}(time(h):time(h+1),:) ), abs(R_wrist_moment{u,1}(time(h):time(h+1),:) ) ];
  le1=size(mom,1);
  le2=size(mom,2);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ANGOLI     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  angles=[abs(Back_angle{u,1}(time(h):time(h+1),1)), abs(L_shoulder_angle{u,1}(time(h):time(h+1),1)), abs(R_shoulder_angle{u,1}(time(h):time(h+1),1)), abs(L_elbow_angle{u,1}(time(h):time(h+1),1)), abs(R_elbow_angle{u,1}(time(h):time(h+1),1)), abs(L_wrist_angle{u,1}(time(h):time(h+1),1)),...
         abs(R_wrist_angle{u,1}(time(h):time(h+1),1)), abs(L_hip_angle{u,1}(time(h):time(h+1),1)), abs(R_hip_angle{u,1}(time(h):time(h+1),1)),abs(L_knee_angle{u,1}(time(h):time(h+1),1)),abs(R_knee_angle{u,1}(time(h):time(h+1),1)),abs(L_ankle_angle{u,1}(time(h):time(h+1),1)),abs(R_ankle_angle{u,1}(time(h):time(h+1),1)) ];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    VELOCITà     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  veloc=[abs(Back_vel{u,1}(time(h):time(h+1),1)), abs(L_shoulder_vel{u,1}(time(h):time(h+1),1)), abs(R_shoulder_vel{u,1}(time(h):time(h+1),1)), abs(L_elbow_vel{u,1}(time(h):time(h+1),1)), abs(R_elbow_vel{u,1}(time(h):time(h+1),1)), abs(L_wrist_vel{u,1}(time(h):time(h+1),1)),...
            abs(R_wrist_vel{u,1}(time(h):time(h+1),1)), abs(L_hip_vel{u,1}(time(h):time(h+1),1)), abs(R_hip_vel{u,1}(time(h):time(h+1),1)),abs(L_knee_vel{u,1}(time(h):time(h+1),1)),abs(R_knee_vel{u,1}(time(h):time(h+1),1)),abs(L_ankle_vel{u,1}(time(h):time(h+1),1)),abs(R_ankle_vel{u,1}(time(h):time(h+1),1)) ];
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    ACCELERAZIONI    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  accel=[abs(Back_acc{u,1}(time(h):time(h+1),1)), abs(L_shoulder_acc{u,1}(time(h):time(h+1),1)), abs(R_shoulder_acc{u,1}(time(h):time(h+1),1)), abs(L_elbow_acc{u,1}(time(h):time(h+1),1)), abs(R_elbow_acc{u,1}(time(h):time(h+1),1)), abs(L_wrist_acc{u,1}(time(h):time(h+1),1)),...
         abs(R_wrist_acc{u,1}(time(h):time(h+1),1)), abs(L_hip_acc{u,1}(time(h):time(h+1),1)), abs(R_hip_acc{u,1}(time(h):time(h+1),1)),abs(L_knee_acc{u,1}(time(h):time(h+1),1)),abs(R_knee_acc{u,1}(time(h):time(h+1),1)),abs(L_ankle_acc{u,1}(time(h):time(h+1),1)),abs(R_ankle_acc{u,1}(time(h):time(h+1),1)) ];
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  arm_fase= abs(ARM{u,1}(time(h):time(h+1),1));
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  momX=zeros(le1,le2/3);
  momnet=zeros(le1,le2/3);
  for j=1:le2/3
   momX(:,j)= mom(:,j*3-2);
   datatot{u,1} { 2+18*(h-1),j+1 } = max(momX(:,j));   %peak_X_mom(u,j)= max(momX(:,j));
   
   momnet(:,j)=   sqrt( mom(:,j*3-2).^2 + mom(:,j*3-1).^2 + mom(:,j*3).^2 );
   datatot{u,1} { 3+18*(h-1),j+1 } = max(momnet(:,j));  %peak_TOT_mom(u,j)= max(momnet(:,j));
  end 
  
  cumx= sum(momX)./100;
  cumtot= sum(momnet)./100;
  for j=1:le2/3
   datatot{u,1} { 4+18*(h-1),j+1 }  = cumx(j);      %cum_X_mom(u,:)= sum(momX)./100;
   datatot{u,1} { 5+18*(h-1), j+1 }  = cumtot(j);   %cum_TOT_mom(u,:)= sum(momnet)./100;
  end    
     
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%valore dei picchi dei momenti negli intervalli considerati 
  pmbackmom= datatot{u,1} { 3+ 18*(h-1), 2 };
  pmLshmom=  datatot{u,1} { 3+ 18*(h-1), 3 };
  pmRshmom=  datatot{u,1} { 3+ 18*(h-1), 4 };
  
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %CALCOLO GLI ARM 
  datatot{u,1} { 18+ 18*(h-1), 15 }= max(arm_fase);
  datatot{u,1} { 19+ 18*(h-1), 15 }= arm_fase(momnet(:,1)==pmbackmom,1);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  for j=1:size(angles,2)
    datatot{u,1} { 6+18*(h-1),j+1 }  =  max(angles(:,j));   %peak_X_ang(u,:)= max(angles); 
    datatot{u,1} { 7+18*(h-1), j+1 }  = max(veloc(:,j));     %peak_X_vel(u,:)= max(veloc);
    datatot{u,1} { 8+18*(h-1), j+1 }  = max(accel(:,j));      %peak_X_acc(u,:)= max(accel);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PM CONSIDERO IL TOT E NON SOLO LA X DEL MOMENTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    datatot{u,1} { 9+ 18*(h-1), j+1 }  = angles(momnet(:,1)==pmbackmom,j);         %BACK         ANGLE
    datatot{u,1} { 10+18*(h-1), j+1 }  = angles(momnet(:,2)==pmLshmom,j);         %LEFT SH
    datatot{u,1} { 11+18*(h-1), j+1 }  = angles(momnet(:,3)==pmRshmom,j);         %RIGHT SH
    
    datatot{u,1} { 12+18*(h-1), j+1 }  = veloc(momnet(:,1)==pmbackmom,j);          %BACK           VEL
    datatot{u,1} { 13+18*(h-1), j+1 }  = veloc(momnet(:,2)==pmLshmom,j);          %LEFT SH
    datatot{u,1} { 14+18*(h-1), j+1 }  = veloc(momnet(:,3)==pmRshmom,j);          %RIGHT SH
    
    datatot{u,1} { 15+18*(h-1), j+1 }  = accel(momnet(:,1)==pmbackmom,j);          %BACK           ACC
    datatot{u,1} { 16+18*(h-1), j+1 }  = accel(momnet(:,2)==pmLshmom,j);          %LEFT SH
    datatot{u,1} { 17+18*(h-1), j+1 }  = accel(momnet(:,3)==pmRshmom,j);          %RIGHT SH
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  end    
         
  
  
   
  end
  
 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADESSO RIEMPIO LA TABELLA DATATOT CON LE ALTRE INFORMAZIONI  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%informazioni specifiche di ogni motion
zz1= Z1{u,1};
zz2= Z2{u,1};
scaf1= 2*(zz1<=0.35) + 3*( (zz1>0.35) && (zz1<=0.65))  +  4*( (zz1>0.65) && (zz1<=0.95))  +...
       5*( (zz1>0.95) && (zz1<=1.25))    +  6*( (zz1>1.25) && (zz1<=1.55))  + 7*( zz1>1.55);
scaf2= 2*(zz2<=0.35) + 3*( (zz2>0.35) && (zz2<=0.65))  +  4*( (zz2>0.65) && (zz2<=0.95))  +...
       5*( (zz2>0.95) && (zz2<=1.25))    +  6*( (zz2>1.25) && (zz2<=1.55))  + 7*( zz2>1.55); 
  
in_bheight= scafftable{scaf1,1};
fin_bheight=scafftable{scaf2,1};

datatot{u,3}=in_bheight/100;
datatot{u,4}=fin_bheight/100;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADESSO RIEMPIO LA TABELLA DATATOT CON i tempi  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%informazioni specifiche di ogni motion
datatot{u,6}= (time(2)-time(1))/100;      %lifting
datatot{u,7}= (time(3)-time(2))/100;      %carrying
datatot{u,8}= (time(4)-time(3))/100;      %lowering
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
 
 
 
  
  
  % end-for dei motion files
 end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%inserisco informazioni comuni a tutti i motion
 datatot{1,5}=name_sub;
 datatot{2,5}='ID';
 datatot{3,5}=id_sub;
 datatot{4,5}='SEX';
 datatot{5,5}=sex_sub;
 datatot{6,5}='HEIGHT';
 datatot{7,5}=HEIGHT{1,1};
 datatot{8,5}='WEIGHT';
 datatot{9,5}=WEIGHT{1,1};
 datatot{10,5}='BOXMASS';
 datatot{11,5}=BOXMASS;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
  
 
 
path=string(subdirec(qq).folder);
pars=strsplit(path,'DATI_MATLAB_INTERPOLATI');
path =pars(1,1)+ "RISULTATI_MATLAB"  +pars(1,end);

name= string(subdirec(qq).name);
pars1=strsplit(name,'INTERP.');
name =pars1(1,1)+ "RIS."  +pars1(1,end);

path= path + "\" + name;
path=char(path);
save(path, 'datatot')

%si può cancellare tutto tranne quello che è fuori dal ciclo qq:qqq (questi ultimi non si devono cancellare )
clearvars -except directory iii nametable scafftable i sub subdirec name_sub parsub id_sub sex_sub qq qqq le0  motion totframes datatot u

 end
 
end