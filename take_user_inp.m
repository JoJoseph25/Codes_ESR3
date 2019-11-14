function [init_struct] = take_user_inp(init_struct)
% function [plot_view,plot_save,mat_save] = take_user_inp()
% 
%  This function asks user to input for various parts of analysis
%  like plots viewing, plot saving, .mat file save. 
% 
%   INPUT:  init_struct - Initialize structure that has various information
%                         (struct)
%     
%   OUTPUT: init_struct - structure that adds following items to it:
%               plot_view - Whether to show figures or close them (boolean)
%               plot_save - Whether to save plot as .jpeg (boolean)
%               mat_save - Whether to save .mat file (boolean)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% Valid Inputs

valid_inp = {'Y','y','N','n'}; % all valid inputs

%% DEFAULT VALUES

prompt = 'Do you want to use default values? (Y/N): '; % Message

check = 0; % check  variable

% While Loop until correct input
while check == 0
    % Try for correct input + check
    try
        inp = input(prompt,'s'); % Input
        inp = validatestring(inp,valid_inp); % Input check
        check = 1;
    catch
        disp("Invalid input you are required to enter from the options above!"); % wrong answer message
    end
end

disp(" ");

if inp=='y' || inp=='Y'
    plot_view = false;
    plot_save = true;
    mat_save = true;
else
    %% PLOT VIEW

    prompt = 'Do you want to view the plots? (Y/N): '; % Message

    check = 0; % check  variable

    % While Loop until correct input
    while check == 0
        % Try for correct input + check
        try
            inp = input(prompt,'s'); % Input
            inp = validatestring(inp,valid_inp); % Input check
            check = 1;
        catch
            disp("Invalid input you are required to enter from the options above!"); % wrong answer message
        end
    end

    disp(" ");

    if inp=='y' || inp=='Y'
        plot_view = true;
    else
        plot_view = false;
    end

    %% PLOT SAVE

    prompt = 'Do you want to save the plots? (Y/N): '; % Message

    check = 0; % check  variable

    % While Loop until correct input
    while check == 0
        % Try for correct input + check
        try
            inp = input(prompt,'s'); % Input
            inp = validatestring(inp,valid_inp); % Input check
            check = 1;

        catch
            disp("Invalid input you are required to enter from the options above!"); % wrong answer message
        end
    end

    disp(" ");

    if inp=='y' || inp=='Y'
        plot_save = true;
    else
        plot_save = false;
    end

    %% MAT SAVE

    prompt = 'Do you want to save .mat file? (Y/N): '; % Message

    check = 0; % check  variable

    % While Loop until correct input
    while check == 0
        % Try for correct input + check
        try
            inp = input(prompt,'s'); % Input
            inp = validatestring(inp,valid_inp); % Input check
            check = 1;
        catch
            disp("Invalid input you are required to enter from the options above!"); % wrong answer message
        end
    end

    disp(" ");

    if inp=='y' || inp=='Y'
        mat_save = true;
    else
        mat_save = false;
    end
end

%% ADD TO INITIALIZE STRUCTURE

init_struct.plot_view=plot_view;
init_struct.plot_save=plot_save;
init_struct.mat_save=mat_save;

end