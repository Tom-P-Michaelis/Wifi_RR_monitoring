%%
% File    : RR_Counter_V5.m    
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
% ┌───────────┐
% │RR Stowood        ├───┬───────┬──────┬─ - - -  
% └───────────┘     │           │          │
%                           │            │         │
%                           ▼            ▼          ▼
%                          P001b        P002b       P003b
%                          │ │         │ │         │ │
%                          │ │         │ └─ ─ ┐  │ └─ ─ ┐
%                          │ │         │            
%                          │ │
%                          │ └-──────┐
%                          │             │
%                          ▼             ▼
%                        Timestamps    Predictions
% ________________________________________________________________________
%
% Improvement Aims from previous versions 
% -----



function RR_Stowood = RR_Counter_Stowood_V10(Sampling_Fs, Subjects)

 Fs = Sampling_Fs;
 
%% Looping over all participants

% Initialising RR_Stowood 

RR_Stowood = {};

for Participant = [1: length(Subjects)] 

%% Reading Subject data 

[~, ~, ~, Airflow] = Data_reader_V5(Subjects(Participant),"Airflow");

%% Removing initial timestamps 

Airflow(:,1) = Airflow(:,1) - Airflow(1,1);

    %% Filtering
    
    Air_window = Airflow(:,2);
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
    
    [Air_window_filt , Filter_config] = IIRFILT(Air_window, Sampling_Fs, Filter_config, false, 1);
    
    % removing padding 
    
    Air_window_filt = Air_window_filt( round(10*Sampling_Fs) : length(Air_window_filt) - round(10*Sampling_Fs) );
    
    Air_window = Air_window( round(10*Sampling_Fs) : length(Air_window)- round(10*Sampling_Fs));
    
    %% Plotting filter effect 
    
    % original data 
    orig = Air_window(400*32:430*32);
    
    tiledlayout(2,3)
    
    nexttile([1,2])
    plot([1:length(orig)]./32, orig./(10^4), 'LineWidth', 1.3)
    axis([0,30, -1,1])
    
    ylabel('Airflow')
    xlabel('Time (s)')
    set(gca,'XMinorTick','on','YMinorTick','on')
%     plotting fft 

nexttile([1,1])
Y = fft(orig);

Y = abs(Y);

Y_Half = (Y(1:floor(length(orig)/2),:));
Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(orig)+1); 
plot(Callibrate, Y_Half)
ylabel('power (a.u.)')
xlabel('Frequency (Hz)')

 axis([0,1.5, 0, max(Y_Half)*1.1])


nexttile([1,2])

    fil = Air_window_filt(400*32:430*32);
    
    
    plot([1:length(fil)]./32, fil./(10^4), 'r', 'LineWidth', 1.3)
    axis([0,30, -1,1])

    ylabel('Airflow')
    xlabel('Time (s)')
    set(gca,'XMinorTick','on','YMinorTick','on')
    
   nexttile([1,1])
Y = fft(fil);

Y = abs(Y);

Y_Half = (Y(1:floor(length(fil)/2),:));
Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(fil)+1); 
    
 plot(Callibrate, Y_Half, 'r')
  axis([0,1.5, 0, max(Y_Half)*1.1])
  ylabel('power (a.u.)')
  xlabel('Frequency (Hz)')

       
%% Peak Counting II

    %% Peak counting II 

% ---------------
% Configuration for Onset/Peak detection
% ---------------

 
 
Onset_config.refractory_period = 1;
Onset_config.learn_period      = 30; % no time to learn as windowed
Onset_config.ssf_window        = 0.65; % should be similar to RR 
Onset_config.NDP               = 5;
Onset_config.max_min_amplitude = prctile(Air_window_filt, 95)*0.4./(10^4); % large to avoid overpicking
% Onset_config.max_min_amplitude = prctile(Air_window_filt, 95)*0.48; % large to avoid overpicking 
Onset_config.min_beat_interval = (1 / 0.002); % lower bound - this has no effect
Onset_config.max_beat_interval = (1 / 0.07); % lower bound on frequency 

% ---------------
% Configuration for Beat SQI
% ---------------

% time vector
TT = [1:length(Air_window_filt)]./32 ;

plot_detail_flag = false;
plot_summary_flag = false;

Beat_onsets = [];

% counting 
[Beat_onsets, Beat_peaks ] = pe.sigproc.pulsew.onset_detect(   ...
        Air_window_filt./(10^4), TT, Fs, Onset_config, plot_detail_flag, plot_summary_flag  ...
        );

    
    %% windowing 
    
 % windowing 
% 30 second overlapping windows 

Overlap = 1 ; 
Win_length = 30;

RR_store = [];
Count2 = 1;

for m = [1: Sampling_Fs*Overlap : length(Airflow)-Sampling_Fs* Win_length ] 
    
    N = 1;
    Temp_Beat = [];
    
    Start = m;
    
    Finish = m+Sampling_Fs* Win_length;
    
    % loop through onsets     
    for q = [1:length(Beat_onsets{1})]
        
        if Beat_onsets{1}(q) >= Start && Beat_onsets{1}(q) <= Finish
        
            % stroing results 
            Temp_Beat(N) = Beat_onsets{1}(q);
            
            N = N + 1 ;
        end    
        
    end 
   
    % checking not empty 
    if length(Temp_Beat) > 1 
    
        % calculating RR 
    Time_wind = (Temp_Beat(length(Temp_Beat)) - Temp_Beat(1))/Fs ;
    Beat_count = length(Temp_Beat) - 1;
    
    RR_pred = 60*Beat_count/Time_wind;
    
    else 
       RR_pred = 0 ;
    end
    
    
 
    %% storing RR and their time stamps 

    RR_store(Count2,:) = [(m+(Win_length*32/2))/32 ,  RR_pred];
    Count2 = Count2 +1 ;
    
end  
    

%% Median Filtering results 

 RR_store =  medfilt1(RR_store , 20,  'truncate') ;

%% Storing Predictions for Participant

RR_Stowood{Participant} = RR_store ; 

end 

end 