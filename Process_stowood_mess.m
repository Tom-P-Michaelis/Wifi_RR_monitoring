%% Importing data 
clc
clear all 

addpath('packages');

%participant number 
Participant = 'P008b';


[~, ~, ~, Airflow] = Data_reader_V5(Participant,"Airflow");

plot([1:length(Airflow(:,2))]./32, Airflow(:,2))

%% processing 

fs = 32 ; 
ppg = Airflow(:,2);
t =  (0:length(Airflow)-1)' ./ fs; 
ppg_t =  (0:length(Airflow)-1)' ./ fs; 
ppg_fs = fs;
ppg_sqi = ones( size(Airflow) );

%% Filter

close all 

Filter_config.method  = {'IIR'; 'IIR'};
Filter_config.type    = {'lowpass'; 'highpass'};
Filter_config.order   = [100, 100];
Filter_config.fc1     = [0.6 , 0.08]; % adjust based on filter results 
Filter_config.fc2     = [0.6, 0.08];

plot_detail_flag = true;

% if filter is defined 
if isempty(Filter_config)
    compute_filter = false;
else
    compute_filter = true;
end

subplot_title = '';

if compute_filter
    ppg = pe.sigproc.filter.filter( ...
        ppg, ppg_fs, Filter_config, t, plot_detail_flag ...
        );
    subplot_title = [subplot_title 'Filt'];
else
    ppg = detrend(ppg);
end

% Invert signal

invert_ppgi = false;

plot_summary_flag = true;
 
if invert_ppgi
    ppg = ppg * -1;
    subplot_title = [subplot_title '+Inv'];
end

fig_row_index = 1;

if plot_summary_flag
    
    fig_row_index = fig_row_index + 1;
    
       [ts_axis_h, fft_axis_h] = ...
    pe.sigproc.pulsew.utils.subplot_signal_in_figure( ...
        figure_h        , ...
        ts_axis_h       , ...
        fft_axis_h      , ...
        num_of_rows     , ...
        num_of_ts_cols  , ...
        num_of_fft_cols , ...
        fig_row_index   , ...
        x_label         , ...
        subplot_title   , ...
        ppg             , ...
        t_plot          , ...
        ppg_fs          , ...
        CP_location     , ...
        Beat_onsets     , ...
        Beat_peaks      , ...
        Beat_sqi      , ...
        plot_valid_data_only , ...
        plot_fft      , ...
        plot_channel  , ...
        sqi_threshold , ...
        y_limit         ...
        );  
    
end

