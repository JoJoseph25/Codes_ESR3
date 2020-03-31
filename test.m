
% function [Exp_Table] = subj_file_save(sub_file,x_num)
% 
%  This function saves subject table and combines all
%  subject struct into one.
% 
%   INPUT:  sub_file - Subject file path (string)
%           x_num - subject ID (string)
%     
%   OUTPUT: Exp_Table - Subject Data combined information (table)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%%  FILENAME CREATION

exp_workdir = 'G:\Socrates\Exp_Data\Henrique_data\Exp_data';
% sub_file = 'G:\Socrates\Exp_Data\Henrique_data\Exp_data\1005';
%
% sub_id = '1005';

subject_num={...
             '1005';'1006';'1007';'1008';...
             '1009';'1010';'1011';'1012';...
             '1013';'1014';'1015';'1016';...
             '1017';'1018';'1019';'1020';...
             '1021';'1022';'1023';'1024';...
             '1025';'1026';'1027';'1028';...
             };
cd(exp_workdir);

%% LOOP OVER ALL SUBJECTS
for i=1:length(subject_num)

%     trial_file = strcat(sub_id,'_',exp_id{i});

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
%          continue
    end
    
    if exist('Exp_Table','var') == 1
        Exp_Table = [Exp_Table; temp_table];
    else
        Exp_Table = temp_table;
    end

end
%% SAVE TO EXCEL
cd(exp_workdir);

writetable(Exp_Table,'Exp_data.xlsx');
