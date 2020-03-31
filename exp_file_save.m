function [Exp_Table] = exp_file_save(exp_workdir,subject_num)
% function Exp_Table = exp_file_save(exp_workdir,subject_num)
% 
%  This function saves subject table and combines all
%  subject struct into one.
% 
%   INPUT:  exp_workdir - Experiment Directory (string)
%           subject_num - Subject ID numbers (cell array)
%     
%   OUTPUT: Exp_Table - Experiment data combined information (table)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%%  FILENAME CREATION

cd(exp_workdir);

disp("Experiment Data Excel: ");

%% LOOP OVER ALL SUBJECTS
for i=1:length(subject_num)

    sub_id = subject_num{i};
    sub_file = strcat(exp_workdir,'\',sub_id);
    
    cd(sub_file);

    filename = strcat(sub_id,'_data.xlsx');

    if exist(filename, 'file') == 2
         % File exists.
         temp_table = readtable(filename);
    else
         % File does not exist.
         disp('File Not Exist')
         continue
    end
    
    % Combine Experiment Data
    if exist('Exp_Table','var') == 1
        Exp_Table = [Exp_Table; temp_table];
    else
        Exp_Table = temp_table;
    end

end

disp('    Experiment Data Created');

%% ASSIGN EMOTION NUM

Exp_Table.emo_num(find(strcmp(Exp_Table.emo_name,'Neutral'))) = 0;
Exp_Table.emo_num(find(strcmp(Exp_Table.emo_name,'Relaxed'))) = 1;
Exp_Table.emo_num(find(strcmp(Exp_Table.emo_name,'Sad'))) = 2;
Exp_Table.emo_num(find(strcmp(Exp_Table.emo_name,'Happy'))) = 3;
Exp_Table.emo_num(find(strcmp(Exp_Table.emo_name,'Fear'))) = 4;

%% SAVE TO EXCEL
cd(exp_workdir);

writetable(Exp_Table,'Exp_data.xlsx');

disp('    Experiment Data Saved');

end