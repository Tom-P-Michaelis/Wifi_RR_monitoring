function [data_filtered, Filter_config] = IIRFILT(data, fs, Filter_config, Adaptor, Time)


 
    
%% Initialising parameters

ppg = data;
t =  (0:length(data)-1)' ./ fs; 
ppg_t =  (0:length(data)-1)' ./ fs; 
ppg_fs = fs;
ppg_sqi = ones( size(data) );

%% Filter parameters 

Filter_config.method  = {'IIR'; 'IIR'};
Filter_config.type    = {'lowpass'; 'highpass'};

if fs == 32 
    
    Filter_config.order   = [8, 8];
    Filter_config.fc1     = [1, 0.08 ];
    Filter_config.fc2     = [1, 0.08 ];

else
    
    Filter_config.order   = [8, 8];
    Filter_config.fc1     = [0.8, 0.08 ];
    Filter_config.fc2     = [0.8, 0.08 ];    
end 

%% Filtering  

% set to true for plot 
if Time >= 300
plot_detail_flag = true;
else 
    plot_detail_flag = false;
end 
plot_detail_flag = false;


% if filter is defined 
if isempty(Filter_config)
    compute_filter = false;
else
    compute_filter = true;
end

subplot_title = '';

if compute_filter
    [ppg, Filter_config] = pe.sigproc.filter.filter( ...
        ppg, ppg_fs, Filter_config, t, plot_detail_flag ...
        );
    subplot_title = [subplot_title 'Filt'];
else
    ppg = detrend(ppg);
end

 

%% Assigning 

data_filtered = ppg;

end 