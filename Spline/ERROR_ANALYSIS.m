%FRAME RANGE ALLARGATO A 7 
%L'ANALISI DEI FRAME DA RIMUOVERE è STATA FATTA USANDO TUTTE LE COMPONENTI
% X-Y-Z DI MOMENTI E ACCELERAZIONE, COSì COME PER L' INTERPOLAZIONE
%considero l'intero esperimento e non taglio alcun frame, ne all'inizio ne
%con l'interpolazione
%INTERPOLAZIONE per momenti, angoli, vel, acc MA SOLO PER: SCHIENA, SPALLE,GOMITI

%------è PERFETTO --------
%AGGIUNGERE IL FOLDER DATI MATLAB E TUTTI I SUBFOLDERS AL PATH
addpath(genpath('D:\SIMONEUTENTE\DESKTOP\PROG_QUA\ESPERIMENTO\QUASISTATIC\DATI_MATLAB'))
clear
close all
clc
directory = dir ('D:\SIMONEUTENTE\DESKTOP\PROG_QUA\ESPERIMENTO\QUASISTATIC\DATI_MATLAB'); % bisogna stare in questa directory 
iii=size(directory,1);
     
 


for vv=3:iii %dentro il FOLDER DATIMATLAB
sub = char( string(directory(vv).folder)+"\"+string(directory(vv).name) );
subdirec = dir (sub);  
qqq=size(subdirec,1);
    
 for pp=3:qqq  %dentro il FOLDER X...X
 load(subdirec(pp).name);%PER OGNI STATIC NELLA CARTELLA
 
 
 %inutili
%  le0=length(Back_acc);
%  motion = cell(le0,1);
%  totframes=zeros(le0,1);
%  for j=1:le0
%   motion{j,1}= ['motion' num2str(j)];
%   totframes(j,1)= size(Back_acc{j,1},1);
%  end
%  datatot=cell(le0,1);

       
       
%comp è il vettore che serve a determinare i range da interpolare,
%NON TOCCARLO
comp= [Back_moment, L_shoulder_moment, R_shoulder_moment, L_elbow_moment, R_elbow_moment,...
       Back_acc,    L_shoulder_acc,     R_shoulder_acc,    L_elbow_acc,   R_elbow_acc];

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% %inizializzo i dati finali  
% ss1= size(Back_moment,1);     ss2= size(Back_moment,2);     z= cell(ss1,ss2);
% %ANGLES  finali
%  Back_angleFIN= z; L_shoulder_angleFIN= z; R_shoulder_angleFIN= z; L_elbow_angleFIN= z; R_elbow_angleFIN= z; 
% %  % L_wrist_angleFIN= z;... R_wrist_angleFIN= z; L_hip_angleFIN= z; R_hip_angleFIN= z; L_knee_angleFIN= z; R_knee_angleFIN= z; L_ankle_angleFIN= z; R_ankle_angleFIN= z; 
% %VELOCITIES  finali
% Back_velFIN= z; L_shoulder_velFIN= z; R_shoulder_velFIN= z; L_elbow_velFIN= z; R_elbow_velFIN= z; 
% %%L_wrist_velFIN= z;...  R_wrist_velFIN= z; L_hip_velFIN= z; R_hip_velFIN= z; L_knee_velFIN= z; R_knee_velFIN= z; L_ankle_velFIN= z; R_ankle_velFIN= z; 
% %ACCELERATIONS  finali
% Back_accFIN= z; L_shoulder_accFIN= z; R_shoulder_accFIN= z; L_elbow_accFIN= z; R_elbow_accFIN= z; 
% % %L_wrist_accFIN= z;... R_wrist_accFIN= z; L_hip_accFIN= z; R_hip_accFIN= z; L_knee_accFIN= z; R_knee_accFIN= z; L_ankle_accFIN= z; R_ankle_accFIN= z; 
% %MOMENTI finali    
%  Back_momentFIN= z;
%  L_shoulder_momentFIN= z;
%  R_shoulder_momentFIN= z;
%  L_elbow_momentFIN= z;
%  R_elbow_momentFIN= z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
clear ANALOG_VIDEO_FRAME_RATIO FILE_NAME  FRAME_RATE Back_moment L_shoulder_moment R_shoulder_moment ...
     L_elbow_moment R_elbow_moment z%BOX_MOVE BOX_PLACE  Back_acc L_shoulder_acc R_shoulder_acc  L_elbow_acc R_elbow_acc

 
 

le=size(comp,1);  

% lista delle coppie di frame da tagliare e interpolare
lista10000=cell(le+1,6);
lista10000(1,2:end)={'Back','L_sh','R_sh','L_el','R_el'};
for j=1:le
  lista10000{j+1,1}= ['motion' num2str(j)];
end

listaDD=lista10000;      %10000 E DD sono 37x6x3 (motion x segmenti x componenti)

listaUNI=lista10000;   %UNI e FIN sono 37x6, solo motion e segmenti su cui interpolare in modo indipendente
listaFIN=lista10000;
% listaMERGE=lista10000;




for u=1:le

%boxtime=[comp(u,11),comp(u,12)];
mom=[comp(u,1),comp(u,2),comp(u,3),comp(u,4),comp(u,5),comp(u,6), comp(u,7),comp(u,8),comp(u,9),comp(u,10)];
nf=size(comp{u,1},1);    %numero frames per ogni motion
nc=size(comp{u,1},2);    %numero di componenti X-Y-Z per ogni dato (è sempre 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TUTTE LE COMPONENTI X-Y-Z PER TUTTA LA DURATA DELL' EXP, 
%non modificare la lunghezza dei frames SE NO GLI EVENTI POI VENGONO
%SBAGLIATI!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
backmom_TOT = cell2mat(mom(1,1));      %backmom_TOT = backmom_TOT(3:nf-2,:);        %backmom_TOT(3:nf-2,:);  % farlo per tom exp, tutta la durata 
lshmom_TOT = cell2mat(mom(1,2));       %lshmom_TOT = lshmom_TOT(3:nf-2,:) ;
rshmom_TOT = cell2mat(mom(1,3));       %rshmom_TOT= rshmom_TOT(3:nf-2,:);
lelbmom_TOT = cell2mat(mom(1,4));      %lelbmom_TOT = lelbmom_TOT(3:nf-2,:) ;
relbmom_TOT = cell2mat(mom(1,5));      %relbmom_TOT= relbmom_TOT(3:nf-2,:);
    
backacc_TOT = cell2mat(mom(1,6));      %backacc_TOT = backacc_TOT(3:nf-2,:);
lshacc_TOT = cell2mat(mom(1,7));       %lshacc_TOT = lshacc_TOT(3:nf-2,:);
rshacc_TOT = cell2mat(mom(1,8));       %rshacc_TOT= rshacc_TOT(3:nf-2,:);
lelbacc_TOT = cell2mat(mom(1,9));      %lelbacc_TOT = lelbacc_TOT(3:nf-2,:);
relbacc_TOT = cell2mat(mom(1,10));     %relbacc_TOT= relbacc_TOT(3:nf-2,:);
    
backdiff_TOT= diff(backmom_TOT);
lshdiff_TOT=  diff(lshmom_TOT);
rshdiff_TOT=  diff(rshmom_TOT);
lelbdiff_TOT= diff(lelbmom_TOT);
relbdiff_TOT= diff(relbmom_TOT);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DATI DA INTERPOLARE ALLA FINE      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INTERPOLATION MOMENT
momGIUNTI= [backmom_TOT, lshmom_TOT, rshmom_TOT, lelbmom_TOT, relbmom_TOT ]; %matrice 1000x 15
momINTERP=momGIUNTI;
%INTERPOLATION ANGLES
angGIUNTI=[Back_angle{u,1}, L_shoulder_angle{u,1}, R_shoulder_angle{u,1}, L_elbow_angle{u,1}, R_elbow_angle{u,1}];    %matrice 1000x 15
%, L_wrist_angle{u,1},...   R_wrist_angle{u,1}, L_hip_angle{u,1}, R_hip_angle{u,1}, L_knee_angle{u,1}, R_knee_angle{u,1}, L_ankle_angle{u,1}, R_ankle_angle{u,1} ];
angINTERP=angGIUNTI;
%INTERPOLATION VELOCITIES
velGIUNTI=[Back_vel{u,1}, L_shoulder_vel{u,1}, R_shoulder_vel{u,1}, L_elbow_vel{u,1}, R_elbow_vel{u,1}];
%, L_wrist_vel{u,1},...     R_wrist_vel{u,1}, L_hip_vel{u,1}, R_hip_vel{u,1}, L_knee_vel{u,1}, R_knee_vel{u,1}, L_ankle_vel{u,1}, R_ankle_vel{u,1} ];%matrice 1000x 15
velINTERP=velGIUNTI;
%INTERPOLATION ACCELERATION
accGIUNTI=[Back_acc{u,1}, L_shoulder_acc{u,1}, R_shoulder_acc{u,1}, L_elbow_acc{u,1}, R_elbow_acc{u,1}];
%, L_wrist_acc{u,1},...    R_wrist_acc{u,1}, L_hip_acc{u,1}, R_hip_acc{u,1}, L_knee_acc{u,1}, R_knee_acc{u,1}, L_ankle_acc{u,1}, R_ankle_acc{u,1} ]; %matrice 1000x 15
accINTERP=accGIUNTI;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



for h=1:3  % per ogni X-Y-Z

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SOLO COMPONENTE ~~~ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
backmom = backmom_TOT (:,h);
lshmom =  lshmom_TOT (:,h);
rshmom =  rshmom_TOT (:,h);
lelbmom = lelbmom_TOT (:,h);
relbmom = relbmom_TOT (:,h);

backacc = backacc_TOT (:,h);
lshacc =  lshacc_TOT (:,h);
rshacc =  rshacc_TOT (:,h);
lelbacc = lelbacc_TOT (:,h);
relbacc = relbacc_TOT (:,h);
  
backdiff= backdiff_TOT(:,h);
lshdiff=  lshdiff_TOT(:,h);
rshdiff=  rshdiff_TOT(:,h);
lelbdiff= lelbdiff_TOT(:,h);
relbdiff= relbdiff_TOT(:,h);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMO CONTROLLO: INDVIDUO FRAMES CON JOINT ACCELERAZIONI OLTRE 10000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
acc=[backacc,lshacc,rshacc,lelbacc,relbacc];

 for j=1:size(acc,2)

  over= find(abs(acc(:,j)) > 10000);   % frames con accel esagerata
  over1 = [1; diff(over)>3; 1];        %distingue il numero di frames consecutivi O CON DISTANZA MINORE DI 3
 
  s=cumsum(over1);              %ai frames consecutivi è assegnato lo stesso numero 
  x =  histc(s,unique(s));       %conto il numero di frames consecutivi
  idx1=x>=3;                  % ind di x di sequenze di frames con almeno 3 elementi, per cui c'è errore solo se 3 almeno frames sono oltre 10000
  idx2=find(x>=3)+1;                % ind di x successivo a idx1
  indicitot = find(over1);          % dentro over trova l'indice del primo frame dopo ogni salto
  indice1=indicitot(idx1);         % dentro over trova l'indice del primo frame di una seq di alm. 3 el 
  indice2=indicitot(idx2)-1;        % dentro over trova l'indice dell'ultimo frame di una seq di alm. 3 el 
  limite1 = over(indice1);          % valore del primo frame della serie 
  limite2 = over( indice2);          %valore dell'ultimo frame della serie  
  lista10000(u+1,j+1,h)=  {[limite1,limite2]};
 end
 clear over over1 s x idx1 idx2 indicitot indice 1 indice2 limite1 limite2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    

 %%                                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SECONDO CONTROLLO: INDVIDUO FRAMES CON DD, il limite è di 5 per tutti
%1 prendo tutti i frames con valore D superiore a abs(5), 
%2 metto nello stesso gruppo i frames con distanza max 10,
%3seleziono solo i gruppi con almeno 2 elementi
% tra questi prendo solo quelli in cui  compaiono D con segni opposti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
der=[backdiff,lshdiff,rshdiff,lelbdiff,relbdiff];

for j=1:size(der,2)
  
  over= find(abs(der(:,j)) > 5);   % frames con accel esagerata
  over1 = [1; diff(over)>10; 1];        %distingue il numero di frames consecutivi O CON DISTANZA MINORE DI 10
  s=cumsum(over1);                  %ai frames consecutivi è assegnato lo stesso numero 
  x =  histc(s,unique(s));          %conto il numero di frames consecutivi
  idx1=x>=2;                  % ind di x di sequenze di frames con almeno 3 elemnti 
  idx2=find(x>=2)+1;                % ind di x successivo a idx1   VEDERE PROT10000 PER OTTIM QUESTA
  indicitot = find(over1);          % dentro over trova l'indice del primo frame dopo ogni salto
  indice1=indicitot(idx1);          % dentro over trova l'indice del primo frame di una seq di alm. 3 el 
  indice2=indicitot(idx2)-1;        % dentro over trova l'indice dell'ultimo frame di una seq di alm. 3 el 
  limite1 = over(indice1);           % valore del primo frame della serie 
  limite2 = over( indice2);          %valore dell'ultimo frame della serie  
  limtot=[limite1,limite2];
   
  ee=zeros(length(limite1),1);
  for i=1:length(ee)
   dd= sign( der(limite1(i):limite2(i),j) );
   ee(i)=max(dd)~=min(dd);
  end
  limtot=limtot(ee==1,:); 
  listaDD(u+1,j+1,h)=  {[limite1,limite2]};
end
 clear over over1 s x idx1 idx2 indicitot indice 1 indice2 limite1 limite2 ee dd
%ADESSO HO I FRAMES CON ALMENO 2 ELEMENTI OPPOSTI CON DISTANZA RELATIVA
%MINORE DI 10
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALISI SPALLA SINISTRA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rangeDD= [  max(1, listaDD{u+1,3,h}(:,1)-5 ),  min( length(acc), listaDD{u+1,3,h}(:,2)+5 )  ];  % range DD in cui trovare errori in un motion file

ee=ones(size(rangeDD,1),1);
for j=1:size(rangeDD,1)     %per ogni intervallo rilevato in un motion file
%adesso se anche un solo frame all'interno di range DD ha un'acc >4000, si
%considera errore
overAsh=  any( abs( lshacc(rangeDD(j,1):rangeDD(j,2)) ) > 4000); 
%adesso se almeno un estremo di tutti quelli oltre 5 per DD gomito sinistro cade in
%tale intervallo, questo è considerato errore.
overDDelb=  any(  rangeDD(j,1)<listaDD{u+1,5,h} &  listaDD{u+1,5,h}<rangeDD(j,2),'all'    ); % se è 0: nessuno nel range o è vuoto

if ~overAsh && ~overDDelb
    ee(j)=0;
%NON SI PUò ALLARGARE IL RANGE ORA, PERCHè LA SPALLA SERVE AL BACK, QUINDI
%SOLO DOPO AVER PROCESSATO IL BACK : listaDD{u+1,4}(j,:)= rangeDD(j,:)
end    
end
listaDD{u+1,3,h}=listaDD{u+1,3,h}(ee==1,:);
clear overAsh overDDelb  ee rangeDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALISI SPALLA DESTRA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rangeDD= [  max(1, listaDD{u+1,4,h}(:,1)-5 ),  min( length(acc), listaDD{u+1,4,h}(:,2)+5 )  ];  % range DD in cui trovare errori in un motion file

ee=ones(size(rangeDD,1),1);
for j=1:size(rangeDD,1)     %per ogni intervallo rilevato in un motion file
%adesso se anche un solo frame all'interno di range DD ha un'acc >4000, si
%considera errore
overAsh=   any( abs( rshacc(rangeDD(j,1):rangeDD(j,2)) ) > 4000);
%adesso se almeno un estremo di tutti quelli oltre 5 per DD gomito dx cade in
%tale intervallo, questo è considerato errore.
overDDelb= any(  rangeDD(j,1)<listaDD{u+1,6,h} &  listaDD{u+1,6,h}<rangeDD(j,2),'all'    );

if ~overAsh && ~overDDelb
    ee(j)=0;
end    
end
listaDD{u+1,4,h}=listaDD{u+1,4,h}(ee==1,:);
clear overAsh overDDelb  ee rangeDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALISI SCHIENA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
rangeDD= [  max(1, listaDD{u+1,2,h}(:,1)-5 ),  min( length(acc), listaDD{u+1,2,h}(:,2)+5 )  ];  % range DD in cui trovare errori in un motion file

ee=ones(size(rangeDD,1),1);
for j=1:size(rangeDD,1)       %per ogni intervallo rilevato in un motion file
%adesso se anche un solo frame all'interno di range DD ha un'acc >5000, si
%considera errore
overAback=   any( abs( backacc(rangeDD(j,1):rangeDD(j,2)) ) > 5000);
%adesso se almeno un estremo di tutti quelli oltre 5 per DD spalle cade in
%tale intervallo, questo è considerato errore.
overDDshl= any(  rangeDD(j,1)<listaDD{u+1,3,h} &  listaDD{u+1,3,h}<rangeDD(j,2),'all'    );
overDDshr= any(  rangeDD(j,1)<listaDD{u+1,4,h} &  listaDD{u+1,4,h}<rangeDD(j,2),'all'    );

if ~overAback && ~overDDshl && ~overDDshr
  ee(j)=0;
end    
end
 listaDD{u+1,2,h}=listaDD{u+1,2,h}(ee==1,:);
 clear overAback overDDshl overDDshr  ee rangeDD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
end  % per ogni X-Y-Z






%%
%FIN QUI LE DUE LISTE 10000 E DD [X] SONO INDIPENDENTI; QUI ORDINO IN COLONNA LE FINESTRE DI FRAMES CON ERRORE E POI FACCIO IL MERGE;
% QUINDI è QUI CHE DEVO AGGIUNGERE LE LISTE 10000 E DD [Y],[Z];  COSì
% 10000[X], 10000[Y],10000[Z], DD[X],DD[Y],DD[Z] SONO TUTTE INDIPENDENTI.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%unisco e ordino le due liste: 10000 e DD
listaDD_X= listaDD(:,:,1);
listaDD_Y= listaDD(:,:,2);
listaDD_Z= listaDD(:,:,3);
lista10000_X= lista10000(:,:,1);
lista10000_Y= lista10000(:,:,2);
lista10000_Z= lista10000(:,:,3);

for j=2:size(listaDD,2)
listaUNI{u+1,j}=  [ listaDD_X{u+1,j}; listaDD_Y{u+1,j}; listaDD_Z{u+1,j} ;...
                    lista10000_X{u+1,j}; lista10000_Y{u+1,j}; lista10000_Z{u+1,j} ]; 
[~,ii]=sort(listaUNI{u+1,j}(:,1));
listaUNI{u+1,j}=listaUNI{u+1,j}(ii,:);
end
clear ii
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%adesso faccio il merge di intervalli vicini che vanno interpolati in un
%solo passaggio
for j=2:size(listaUNI,2)
    for h=size(listaUNI{u+1,j},1):-1:2
       if listaUNI{u+1,j}(h,1) >=  listaUNI{u+1,j}(h-1,1)  &&  listaUNI{u+1,j}(h,1) <=  listaUNI{u+1,j}(h-1,2)+17 
       listaUNI{u+1,j}(h-1,:)=  [listaUNI{u+1,j}(h-1,1),    max(  listaUNI{u+1,j}(h,2), listaUNI{u+1,j}(h-1,2)  )  ];
       listaUNI{u+1,j}(h,:)=[];
       end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%REGOLA: ALLARGO 5--> DIST 13;  ALL.6-->15; ALL.7-->17 


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%adesso allargo ogni finestra di 10 frames per fare un'interpolazione congrua 
for j=2:size(listaUNI,2)
listaFIN{u+1,j}= [  max(1, listaUNI{u+1,j}(:,1)-7 ),  min( length(acc), listaUNI{u+1,j}(:,2)+7 )  ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%listaFIN è LA LISTA DEFINITIVA CONTENETE I GLI INTERVALLI COMPLETI DI
%FRAMES DA INTERPOLARE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 






%%
%INTERPOLAZIONE
%per ogni motion (U), per ogni giunto del corpo umano (J) si esegue la spline interpolante negli intervalli di interesse ([0,1;2,3]--> I),
%bisogna interpolare in tutte e tre le direzioni X-Y-Z (H)

for j=2:size(listaFIN,2)  %per ogni segmento del corpo
    nint= size(listaFIN{u+1,j},1);
    
    
    for i=1:nint            %per ogni intervallo di interesse [......] relativo al seg nel motion, lascio 3 frame prima e dopo come poli fissi
        %I FRAME NON CAMBIANO PER ang,vel,acc; MENTRE I VALORI SI
        %x=frame-polifissi, STESSI per ang-vel..----  xx= x_spline, STESSI
        %qq= valore di x, CAMBIA                ----  yy= valore di xx, CAMBIA 
        x = [  max(1, listaFIN{u+1,j}(i,1)-3 ): listaFIN{u+1,j}(i,1) ,  listaFIN{u+1,j}(i,2) : min( length(acc), listaFIN{u+1,j}(i,2)+3 )   ];
        xx= x(1):x(end);
        yy=   zeros(length(xx), nc);
        yyang=zeros(length(xx), nc);
        yyvel=zeros(length(xx), nc);
        yyacc=zeros(length(xx), nc);
        %%%%%%%%%%%%%%%%%%
        ngiu= (j-1)*3-2;                %indice di momGIUNTI
        qq=    momGIUNTI(x, ngiu:ngiu+2 ); %valore dei momenti nei frame del vettore x
        qqang= angGIUNTI(x, ngiu:ngiu+2 );
        qqvel= velGIUNTI(x, ngiu:ngiu+2 );
        qqacc= accGIUNTI(x, ngiu:ngiu+2 );
        %%%%%%%%%%%%%%%%%
        for h = 1:nc                   %per ogni coordinata del giunto umano
            yy(:,h)    = spline(x,qq(:,h),xx); 
            yyang(:,h) = spline(x,qqang(:,h),xx);
            yyvel(:,h) = spline(x,qqvel(:,h),xx);
            yyacc(:,h) = spline(x,qqacc(:,h),xx);
        end
        
        momINTERP(xx,ngiu:ngiu+2)= yy;
        angINTERP(xx,ngiu:ngiu+2)= yyang;
        velINTERP(xx,ngiu:ngiu+2)= yyvel;
        accINTERP(xx,ngiu:ngiu+2)= yyacc;
    end
    
    Back_moment{u,1} =      momINTERP(:,1:3);
    L_shoulder_moment{u,1}= momINTERP(:,4:6);
    R_shoulder_moment{u,1}= momINTERP(:,7:9);
    L_elbow_moment{u,1}=    momINTERP(:,10:12);
    R_elbow_moment{u,1}=    momINTERP(:,13:15);
   %ANGLES  finali
   Back_angle {u,1}=         angINTERP(:,1:3);
   L_shoulder_angle {u,1}=   angINTERP(:,4:6);
   R_shoulder_angle {u,1}=   angINTERP(:,7:9);
   L_elbow_angle {u,1}=      angINTERP(:,10:12);
   R_elbow_angle {u,1}=      angINTERP(:,13:15);
   %VELOCITIES  finali
   Back_vel {u,1}=          velINTERP(:,1:3);
   L_shoulder_vel {u,1}=    velINTERP(:,4:6);
   R_shoulder_vel {u,1}=    velINTERP(:,7:9);
   L_elbow_vel {u,1}=       velINTERP(:,10:12);
   R_elbow_vel {u,1}=       velINTERP(:,13:15);
   %ACCELERATIONS  finali
   Back_acc {u,1}=        accINTERP(:,1:3);
   L_shoulder_acc {u,1}=  accINTERP(:,4:6);
   R_shoulder_acc {u,1}=  accINTERP(:,7:9);
   L_elbow_acc {u,1}=     accINTERP(:,10:12);
   R_elbow_acc {u,1}=     accINTERP(:,13:15);

    
    
end






% end-for dei motion files
end


path=string(subdirec(pp).folder);
pars=strsplit(path,'DATI_MATLAB');
path =pars(1,1)+ "DATI_MATLAB_INTERPOLATI"  +pars(1,end);

name= string(subdirec(pp).name);
pars1=strsplit(name,'.');
name =pars1(1,1)+ "INTERP."  +pars1(1,end);




path= path + "\" + name;
path=char(path);
save(path, 'Back_moment','L_shoulder_moment','R_shoulder_moment','L_elbow_moment','R_elbow_moment','L_wrist_moment','R_wrist_moment',... 
 'Back_angle','L_shoulder_angle','R_shoulder_angle','L_elbow_angle','R_elbow_angle','L_wrist_angle','R_wrist_angle',...
 'L_hip_angle','R_hip_angle','L_knee_angle','R_knee_angle','L_ankle_angle','R_ankle_angle',...
 'Back_vel','L_shoulder_vel','R_shoulder_vel','L_elbow_vel','R_elbow_vel','L_wrist_vel','R_wrist_vel',...
 'L_hip_vel','R_hip_vel','L_knee_vel','R_knee_vel','L_ankle_vel','R_ankle_vel',...
  'Back_acc','L_shoulder_acc','R_shoulder_acc','L_elbow_acc','R_elbow_acc','L_wrist_acc','R_wrist_acc',...
 'L_hip_acc','R_hip_acc','L_knee_acc','R_knee_acc','L_ankle_acc','R_ankle_acc',...
 'CROSS_LAB_0','CROSS_LAB_1','CROSS_LAB_2','START_LIFTING','END_LIFTING','START_LOWERING','END_LOWERING','BOX_MOVE','BOX_PLACE',...
 'HEIGHT','WEIGHT','BOXMASS','Y1','Y2','Z1','Z2','ARM') 


clearvars  -except directory iii i sub subdirec pp qqq le0  motion totframes datatot u


end
 
end




