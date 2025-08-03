% Example on how to process each subcarrier for each PCA or ICA component,
% individually.
%
% The parameters in this script have been tuned to estimate heart rate,
% you will need to adjust them for respiratory rate.
%
% Follow the main function "pe.sigproc.pulsew.preprocess_waveform" to
% understand the meaning of the parameters and read the referneced
% paper

%
% Run the startup.m script before using these packages


% ---------------
% FILL THE FOLLOWING VARIABLES ACCORDINGLY TO THE RECORD
% YOU ARE PROCESSING
% ---------------

record_name = 'S01a';  % The name you assigned to the record
num_of_frames = 1;     % The n umbrer of sample in the array
colour_channels = 1;   % One single channel per subcarrier and ICA component
fs = 50;               % Sampling frequency

blue_colour        = [0 .4 .7];

% Create the record structure for the channel and ICA componentn

WiFi_record.record_name     = record_name;
WiFi_record.name            = 'Wifi amplitude';
WiFi_record.type            = 'PPG';
WiFi_record.units           = 'a.d.u.';
WiFi_record.source          = 'Rasberry Pi Nexmon';
WiFi_record.video_path      = 'path to the file';
WiFi_record.ts              = zeros(num_of_frames, colour_channels, 'double'); % "ts" means time_series, the raw wifi data
WiFi_record.sqi             = ones( size(WiFi_record.ts) );     % "sqi" is the is the signal quality indices for every single sample in the "ts" array
WiFi_record.fs              = fs;
WiFi_record.t               = (0:num_of_frames-1)' ./ fs;       % "t" is the time vector
WiFi_record.start_frame     = 1;     % You don't need to use this.
WiFi_record.channel_num     = colour_channels;   % you have a single channel
WiFi_record.colour          = blue_colour;       % The colour map that will be used to plot the time series
WiFi_record.channel_labels  = {'subcarrier n'};  % name to show in the plot legend


% ---------------
% Plotting properties and feedback
% ---------------


% 
% Adjust the following variables depending on how many detailed figures 
% you want
plot_summary_flag             = true;
plot_record_summary_flag      = false;
plot_record_detail_flag       = false;
plot_record_extra_detail_flag = false;
plot_filter_design_flag       = false;
plot_channel                  = 1; % Use the green channel for example plots

% Write extract log messages to the console
verbose_flag = true;


% ---------------
% Parallel processing
% ---------------


% Cofigure how many records you want to run in parallel.
% If you enable threading, all the plot detail flags will be disabled
% 
num_of_threads = 0;


%% Image and signal processing paramaters
%
% Adjust these parameters only after you read the papers the functions in
% the "pe.sigproc" package reference
% ========================================================================


% ---------------
% Blobal configuration parameters
% ---------------

PPGi_config.target_fs          = 0;
PPGi_config.invert_ppgi        = true;
PPGi_config.valid_low_freq     = 0.5;  % This the lowest valid frequency for heart rate ranges, modify for RR
PPGi_config.valid_high_freq    = 2.0;  % This the highest valid frequency for heart rate ranges, modify for RR

% ---------------
% Configuration for Change Point detection
% ---------------

CP_config.wlen         = [ 30   15  10 ];
CP_config.window_delta = [ 5     5   5 ];
CP_config.threshold    = 0.8;

% ---------------
% Configuration for Onset/Peak detection
% ---------------

Onset_config.refractory_period = 0.34;
Onset_config.learn_period      = 60;
Onset_config.ssf_window        = 0.17;
Onset_config.NDP               = 2.5;
Onset_config.max_min_amplitude = 0.4;
Onset_config.min_beat_interval = (1 / PPGi_config.valid_high_freq);
Onset_config.max_beat_interval = (1 / PPGi_config.valid_low_freq);

% ---------------
% Configuration for Beat SQI
% ---------------

SQI_config.wlen            = 20; 
SQI_config.window_delta    = 15; 
SQI_config.beat_boundary   = 2;
SQI_config.valid_low_freq  = PPGi_config.valid_low_freq;
SQI_config.valid_high_freq = PPGi_config.valid_high_freq;
SQI_config.template.beat_threshold  = 0.5;
SQI_config.dtw.pla_error_thresholds = [ 1  0.5 0.25 0];
SQI_config.dtw.clip_amp_threshold   = 0.01;
SQI_config.dtw.min_dist_threshold   = 15; %10

% ---------------
% Configuration to compute metrics for the PPGi signal such as
% pulse amplitude, width, etc.
% ---------------

PPGi_metrics_config.compute_metrics_flag = true;

% ---------------
% Configuration for the rough HR estimates
% ---------------

Estimation_config.algorithm_name = 'BEAT_PEAKS';
Estimation_config.min_error_threshold = 1.0;
Estimation_config.max_error_threshold = 10.0;

% ---------------
% Activity/Motion quality index
% ---------------

activity_index = [];

% ---------------
% Configuration for the filtering step
% ---------------

Filter_config.method  = {'FIR'; 'FIR'};
Filter_config.type    = {'highpass'; 'lowpass'};
Filter_config.order   = [250; 150];
Filter_config.fc1     = [0.7, 2.5];
Filter_config.fc2     = [0.7, 2.5];

% Plot the filter design


    
    filter_count = length( Filter_config.type );

    for i = 1 : filter_count
        
        filter_breaks    = [Filter_config.fc1(i) Filter_config.fc2(i) ];
        filter_algorithm = [];
        
        Filter_config.filter_object{i} = pe.sigproc.filter.create_filter(  ...
            filter_breaks               , ...
            fs , ...
            Filter_config.order(i)      , ...
            Filter_config.type(i)       , ...
            Filter_config.method(i)     , ...
            filter_algorithm            , ...
            plot_filter_design_flag       ...
            );
    end



%% Process data


    [ ...
        WiFi_record.ts             , ...
        WiFi_record.t              , ...
        WiFi_record.fs             , ...
        WiFi_record.sqi            , ...
        WiFi_record.CP_location    , ...
        WiFi_record.Beat_onsets    , ...
        WiFi_record.Beat_peaks     , ...
        WiFi_record.Beat_sqi       , ...
        WiFi_record.Beat_t           ...
    ] = ...
pe.sigproc.pulsew.preprocess_waveform( ...
        WiFi_record.ts                , ...
        WiFi_record.t                 , ...
        WiFi_record.fs                , ...
        PPGi_config.target_fs      , ...
        PPGi_config.invert_ppgi    , ...
        Filter_config              , ...
        CP_config                  , ...
        Onset_config               , ...
        SQI_config                 , ...
        activity_index             , ...
        plot_detail_flag           , ...
        plot_extra_detail_flag     , ...
        plot_channel                 ...
    );

