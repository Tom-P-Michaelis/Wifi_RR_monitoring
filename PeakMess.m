%% Importing data 

clear all 
close all 

%participant number 
Participant = 'P004b';

[~, ~, ~, Airflow] = Data_reader_V5(Participant,"Airflow");

addpath('C:\Users\tommi\OneDrive\Documents\4th Year Eng\4YP\Matlab\packages')

 
%% sudo windowing 

% Stowood sampling frequency 
Fs = 32 ;
Sampling_Fs = Fs ;


Airflow3 = (Airflow(60*Fs:90*Fs,2));
Airflow4 = (Airflow(150*Fs:180*Fs,2));
Airflow1 = (Airflow(1150*Fs:1180*Fs,2));
Airflow2 = (Airflow(100*Fs:130*Fs,2));



%% Filtering 

%%  Filtering Adaptively 

for Air_window = [ Airflow1, Airflow3, Airflow4, Airflow2 ]
    
    %% detrend
    Air_window = Air_window - mean(Air_window);
    Air_window = detrend(Air_window);
    
    %% Filtering
    
    % initializing filter config 
    Filter_config.a = 1 ;
   
    % finding start value 
    Mem1 = median(Air_window(1:1*round(Sampling_Fs)));
    
    % finding end value 
    Mem2 = median( Air_window(length(Air_window)-1*round(Sampling_Fs):length(Air_window)) );
    
    % creating padding array that is 10 seconds long 
    
    Mem1 = Mem1 * ones(round(10*Sampling_Fs),1)  ;
    Mem2 = Mem2 * ones(round(10*Sampling_Fs),1)  ;
    
    % adding padding array 
    Air_window = [Mem1; Air_window; Mem2];
    
    [Air_window_filt , Filter_config] = IIRFILT(Air_window, Sampling_Fs, Filter_config, false);
    
    % removing padding 
    
    Air_window_filt = Air_window_filt( round(10*Sampling_Fs) : length(Air_window_filt) - round(10*Sampling_Fs) );
    
    Air_window = Air_window( round(10*Sampling_Fs) : length(Air_window)- round(10*Sampling_Fs));
    
    %% plotting 
    
    % plot airflow, airflow filt 
    figure
    plot([1:length(Air_window_filt)], Air_window_filt, [1:length(Air_window)], Air_window) 
    legend('filtered', 'non - filtered' )
    
    
    %% Peak Counting 
    
    %% Peak counting II 

% ---------------
% Configuration for Onset/Peak detection
% ---------------
 
Onset_config.refractory_period = 0.05;
Onset_config.learn_period      = 0; % no time to learn as windowed
Onset_config.ssf_window        = 0.3; % should be similar to RR 
Onset_config.NDP               = 4;
Onset_config.max_min_amplitude = prctile(Air_window_filt, 95)*0.45; % large to avoid overpicking 
Onset_config.min_beat_interval = (1 / 0.002); % lower bound - this has no effect
Onset_config.max_beat_interval = (1 / 0.2); % upper bound 

% ---------------
% Configuration for Beat SQI
% ---------------

TT = [1:length(Air_window_filt)]./32 ;

plot_detail_flag = true;

Beat_onsets = [];

% counting 
[Beat_onsets, Beat_peaks ] = pe.sigproc.pulsew.onset_detect(   ...
        Air_window_filt, TT, Fs, Onset_config, plot_detail_flag, true  ...
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
pause

    
    
    
    
end 
 

 %% 
 length(Beat_onsets{1})





