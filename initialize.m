function [init_struct] = initialize()
% function [init_struct] = initialize()
% 
%  This function declares sampling frequency and filter coefficients
%  used to clean the data from force plate and marker data. 
%    
%   INPUT:  NA
%     
%   OUTPUT: init_struct - structure that contains following items:
%               Fc - Frequncy of recording (scalar)
%               dt - time difference between consequetive data points (scalar) 
%               b_fp - numerator filter coefficient for forceplate data (scalar)
%               a_fp - denominator filter coefficient for forceplate data (scalar)
%               b_mrk - numerator filter coefficient for marker data (scalar)
%               a_mrk - denominator filter coefficient for marker data (scalar)
%             
% written by Joel V Joseph (josephjo@post.bgu.ac.il)

%% SAMPLING FREQUENCY

Fs =120; % 120 Hz for both marker and force plate data
dt=1/Fs; % 120Hz

%% FORCEPLATE DATA

Fc =20; % Cutoff frequency for force plate 20Hz
order=2; % Order for buttworth filter

%Buttworth filter design to get filter coefficients
[b_fp,a_fp] = butter(order,2*Fc/Fs);

%% MARKER DATA

Fc_mrk =10; % Cutoff frequency for force plate 20Hz
order=2; % Order for buttworth filter

%Buttworth filter design to get filter coefficients
[b_mrk,a_mrk] = butter(order,2*Fc_mrk/Fs);

%% INITIALIZE STRUCTURE

init_struct=struct; % initialization structure

init_struct.Fs=Fs;
init_struct.dt=dt;
init_struct.b_fp=b_fp;
init_struct.a_fp=a_fp;
init_struct.b_mrk=b_mrk;
init_struct.a_mrk=a_mrk;

%% DISPLAY

disp('Initalization Done'); % print out to command line

end