function [Subj_Table] = subj_file_save(sub_file,x_num)
% function [Subj_Table] = subj_file_save(sub_file,x_num)
% 
%  This function saves subject table and combines all
%  subject struct into one.
% 
%   INPUT:  sub_file - Subject file path (string)
%           x_num - subject ID (string)
%     
%   OUTPUT: Subj_Table - Subject Data combined information (table)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%%  FILENAME CREATION

sub_id = x_num;
exp_id = {'dy4';'dy5';'dy6';'dy7';'dy8'};

cd(sub_file);

disp("Subject Data Excel: ");

%% QUESTIONNARE DATA
sub_num = str2num(sub_id(3:4));
emo_workdir = '../../Questionnaire';

filename_question = 'motion_emotion.csv';

addpath(emo_workdir);

opts_question = detectImportOptions(filename_question,'FileType','text');
opts_question.VariableNamesLine=1; % Where the variable names are located
opts_question.DataLine=4; % Line number where the table is located

Question_Data = readtable(filename_question,opts_question);
sub_row = find(Question_Data.SN == sub_num);

age = Question_Data.Q69(sub_row);
gender = Question_Data.Q70(sub_row);

ER_1 = Question_Data.Q49_1(sub_row);
ER_2 = Question_Data.Q51_1(sub_row);
ER_3 = Question_Data.Q53_1(sub_row);
ER_4 = Question_Data.Q55_1(sub_row);
ER_5 = Question_Data.Q57_1(sub_row);
ER_6 = Question_Data.Q59_1(sub_row);
ER_7 = Question_Data.Q61_1(sub_row);
ER_8 = Question_Data.Q63_1(sub_row);
ER_9 = Question_Data.Q65_1(sub_row);
ER_10 = Question_Data.Q67_1(sub_row);

ER_cog_reapp = mean([ER_1 ER_3 ER_5 ER_7 ER_8 ER_10]);
ER_emo_sup = mean([ER_2 ER_4 ER_6 ER_9]);

RF_1 = Question_Data.Q42_1(sub_row);
RF_2 = Question_Data.Q42_2(sub_row);
RF_3 = Question_Data.Q42_3(sub_row);
RF_4 = Question_Data.Q42_4(sub_row);
RF_5 = Question_Data.Q42_5(sub_row);
RF_6 = Question_Data.Q42_6(sub_row);
RF_7 = Question_Data.Q42_7(sub_row);
RF_8 = Question_Data.Q42_8(sub_row);
RF_9 = Question_Data.Q42_9(sub_row);
RF_10 = Question_Data.Q42_10(sub_row);
RF_11 = Question_Data.Q42_11(sub_row);
RF_12 = Question_Data.Q42_12(sub_row);
RF_13 = Question_Data.Q42_13(sub_row);
RF_14 = Question_Data.Q42_14(sub_row);

RF_promotion = mean([RF_1 RF_2 RF_3 RF_4 RF_5 RF_6 RF_7]);
RF_prevention = mean([RF_8 RF_9 RF_10 RF_11 RF_12 RF_13 RF_14]);

rmpath(emo_workdir);

disp('    Questionnare Data Fetched');

%% LOOP OVER ALL TRIALS
for i=1:length(exp_id)

    trial_file = strcat(sub_id,'_',exp_id{i});

    cd(strcat(sub_file,'\',trial_file));

    filename = strcat(trial_file,'_struct.xlsx');
    
    % Check if file exist
    if exist(filename, 'file') == 2
         % File exists.
         temp_table = readtable(filename);
    else
         % File does not exist.
         disp('    File Not Exist')
         continue
    end
    
    % Manipulation Check
    if strcmp(temp_table.emo_name,'Neutral')
        familiarity =0 ;
        manpul_sad = 0;
        manpul_relax = 0;
        manpul_happy = 0;
        manpul_fear = 0;
        valence = 0;
        arousal = 0;
    elseif strcmp(temp_table.emo_name,'Relaxed')
        familiarity = Question_Data.Q97_1(sub_row); 
        manpul_sad = Question_Data.Q98_1(sub_row);
        manpul_relax = Question_Data.Q98_2(sub_row);
        manpul_happy = Question_Data.Q98_3(sub_row);
        manpul_fear = Question_Data.Q98_4(sub_row);
        valence = Question_Data.Q75_1(sub_row);
        arousal = Question_Data.Q77(sub_row);
    elseif strcmp(temp_table.emo_name,'Sad')
        familiarity = Question_Data.Q99(sub_row); 
        manpul_sad = Question_Data.Q100_3(sub_row);
        manpul_relax = Question_Data.Q100_4(sub_row);
        manpul_happy = Question_Data.Q100_5(sub_row);
        manpul_fear = Question_Data.Q100_6(sub_row);
        valence = Question_Data.Q78(sub_row);
        arousal = Question_Data.Q79_1(sub_row);
    elseif strcmp(temp_table.emo_name,'Happy')
        familiarity = Question_Data.Q101(sub_row); 
        manpul_sad = Question_Data.Q102_3(sub_row);
        manpul_relax = Question_Data.Q102_4(sub_row);
        manpul_happy = Question_Data.Q102_5(sub_row);
        manpul_fear = Question_Data.Q102_6(sub_row);
        valence = Question_Data.Q80(sub_row);
        arousal = Question_Data.Q81(sub_row);
   elseif strcmp(temp_table.emo_name,'Fear')
        familiarity = Question_Data.Q103(sub_row); 
        manpul_sad = Question_Data.Q104_3(sub_row);
        manpul_relax = Question_Data.Q104_4(sub_row);
        manpul_happy = Question_Data.Q104_5(sub_row);
        manpul_fear = Question_Data.Q104_6(sub_row);
        valence = Question_Data.Q82(sub_row);
        arousal = Question_Data.Q83(sub_row);
    end
    
    % Create and concatenate table
    if exist('Subj_Table','var') == 1
        append_table = [array2table([age gender familiarity manpul_sad manpul_relax manpul_happy manpul_fear valence arousal ER_cog_reapp ER_emo_sup RF_promotion RF_prevention]) temp_table];
        append_table.Properties.VariableNames([1 2 3 4 5 6 7 8 9 10 11 12 13]) = {'Age','Gender','Familiarity','Sad','Relax','Happy','Fear','Valence','Arousal','ER_Cog_Reapp','ER_Emo_Supp','RF_promote','RF_prevent'};
        Subj_Table = [Subj_Table; append_table];
    else
        Subj_Table = [array2table([age gender familiarity manpul_sad manpul_relax manpul_happy manpul_fear valence arousal ER_cog_reapp ER_emo_sup RF_promotion RF_prevention]) temp_table];
        Subj_Table.Properties.VariableNames([1 2 3 4 5 6 7 8 9 10 11 12 13]) = {'Age','Gender','Familiarity','Sad','Relax','Happy','Fear','Valence','Arousal','ER_Cog_Reapp','ER_Emo_Supp','RF_promote','RF_prevent'};
    end

end

disp('    Subject Data Created');

%% SAVE TO EXCEL
cd(sub_file);

Subj_Table = [Subj_Table(:,14:17)  Subj_Table(:,1:13) Subj_Table(:,18:end)];

writetable(Subj_Table, strcat(sub_id,'_data.xlsx'));

disp('    Subject Data Saved');

end