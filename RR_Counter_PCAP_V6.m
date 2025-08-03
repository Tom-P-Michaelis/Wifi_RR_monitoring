%%
% File    : RR_Counter_PCAP_V5.m    
% Author  : Tom Michaelis (tom.michaelis.tm@googlemail.com)
% Created : 09/04/2021
% Updated : 
% ________________________________________________________________________
% 
% This file is part of Wi-Fi Vitals Monitoring
%
% Wi-Fi Vitals Monitoring: A library of MATLAB scripts for predicting 
% physiological measurements from CSI Wi-Fi data. 
%
% ________________________________________________________________________
%
% DESCRIPTON:
% ----------
%
%  RR_Counter_V5.m is a simple peak counting algorithm that predicts the
%  respiratory rate from Stowood nasal pressure canula data. 
%
% Reference:
%
% (https://link.springer.com/content/pdf/10.1007/s10439-007-9428-1.pdf)
%
% INPUT
% -----
%
% Subjects: The list of subject numbers to be included in the analysis 
% 
% Fs: The recoding frequency of the Stowood device. Currently 32Hz. 
% 
% OUTPUT
% -----
% 
% RR_Stowood: A cell vector containing the results of the RR predictions:
%
% ________________________________________________________________________
%
% Improvement Aims from previous versions 
% -----
% 
% > Improve peak count to detect failure 
% > Combine for all subjects [DONE]
% > Change to power based filter selection
% > Introduce results flagging 

function Predictions = RR_Counter_PCAP_V6(Airflow, Fs, Time)

    
    Air_window_filt = Airflow;
   

    %% Peak counting II 

% ---------------
% Configuration for Onset/Peak detection
% ---------------
 
Onset_config.refractory_period = 0.6;
Onset_config.learn_period      = 10; % no time to learn as windowed
Onset_config.ssf_window        = 0.6; % should be similar to RR 
Onset_config.NDP               = 4;
Onset_config.max_min_amplitude = prctile(Air_window_filt, 95)*0.48; % large to avoid overpicking 
Onset_config.min_beat_interval = (1 / 0.002); % lower bound - this has no effect
Onset_config.max_beat_interval = (1 / 0.2); % upper bound 

% ---------------
% Configuration for Beat SQI
% ---------------

TT = [1:length(Air_window_filt)]./32 ;

plot_detail_flag = false;
plot_summary_flag = false;

% Plotting commands
% if Time > 80
%     plot_detail_flag = true;
%     plot_summary_flag = true;
% end

Beat_onsets = [];

% counting 
[Beat_onsets, Beat_peaks ] = pe.sigproc.pulsew.onset_detect(   ...
        Air_window_filt, TT, Fs, Onset_config, plot_detail_flag, plot_summary_flag  ...
        );

    % checking not empty 
    if length(Beat_onsets{1}) > 1 
    
        % calculating RR 
    Time_wind = (Beat_onsets{1}(length(Beat_onsets{1}),1) - Beat_onsets{1}(1,1))/Fs ;
    Beat_count = length(Beat_onsets{1}) - 1;
    
    RR_pred = 60*Beat_count/Time_wind;
    
    else 
       RR_pred = 0 ;
    end 
    
%     if Time > 80
% pause
% end
    
   
    
%% storing RR and their time stamps 


Predictions = RR_pred;
end 