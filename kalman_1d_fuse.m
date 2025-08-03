%%
% File    : kalman_1d.m    
% Author  : Mauricio Villarroel
% Created : Feb 13 2014
% Updated : $Id$
% ________________________________________________________________________
%
% This file is part of ParamEstimator
%
% ParamEstimator: Library for deriving physiological parameters from cameras
% Copyright (C) 2011 Mauricio Villarroel. All rights reserved.
%
% PROPRIETARY/CONFIDENTIAL. Use is subject to license terms.
%
% You may contact the author by e-mail (villarroel DOT mauricio AT gmail)
% or postal mail to:
%
%    - Mauricio Villarroel
%    - Room 20.62
%    - Institute of Biomedical Engineering
%    - Department of Engineering Science
%    - University of Oxford
%    - Old Road Campus Research Building (off Roosevelt Drive),
%    - Oxford OX3 7DQ
%    - United Kingdom
%
% ________________________________________________________________________
%
% DESCRIPTON:
% ----------
%
%    Simple Kalman filter for 1D time series. From the referenced papers,
%  it assumes that:
%
%  A = 1
%  H = 1
%  B = 0
%  u = 0
%
% Reference:
%
% * Welch, Greg, and Gary Bishop. "An introduction to the Kalman filter."
%  (1995).
%
% * Li, Qiao, Roger G. Mark, and Gari D. Clifford. 
%   "Robust heart rate estimation from multiple asynchronous noisy 
%   sources using signal quality indices and a Kalman filter." 
%   Physiological measurement 29.1 (2008): 15.
%
% INPUT
% -----
%
%    Z      : 1D Column-based input time series
%
%   sqi     : SQI for the input data
%
%   fs      : Sample frequency
%
%   min/max_error_threshold  : Minimum and maximum allowed error in the
%             a posteriori error covariante "P" before the
%             SQI is assigned 0
%
% OUTPUT
% -----
%
% ________________________________________________________________________

function [ ...
            X      , ...
            X_sqi  , ...
            r      , ...
            P        ...
         ] = ...
    kalman_1d_fuse( ...
            Zs                   , ...
            sqis                 , ...
            fs                  , ...
            min_error_threshold , ...
            max_error_threshold , ...
            plot_summary_flag   , ...
            plot_detail_flag    , ...
            t                     ...
        )

narginchk(5, inf);

%% Process function parameters with default values
%

if nargin < 6
    plot_summary_flag = false;
end
if nargin < 7
    plot_detail_flag = false;
end
if nargin < 8
    t = [];
end

%% Initializing the inputs as candidate 1

Z = Zs{1}(:,2);
sqi = sqis{1};


%% Global configurable parameters
%

[signal_len, channels] = size(Z);

learning_period  = 60;
learn_period_len = min( ceil(learning_period * fs), signal_len );
num_of_learning_segments = 10;

sqi_threshold = 0.8;

Xa = zeros( signal_len + 1, 1 );  % Apriori state estimates
X  = zeros( signal_len,     1 );  % Aposteriori state estimates
Pa = zeros( signal_len + 1, 1 );  % Apriori estimate error covariance 
P  = zeros( signal_len,     1 );  % Aposteriori estimate error covariance
K  = zeros( signal_len,     1 );  % Kalman gain

Q  = ones( signal_len, 1 );   % Process noise variance estimates
                              % Bishop, Eq. 1.3
                                     
R  = ones( signal_len, 1 );   % Measurement noise variance estimates
                              % Bishop, Eq. 1.4

r  = ones( signal_len, 1 );   % Kalman innovation/residual

%% Checks

if channels ~= 1 
    error( 'Input signal MUST have one data channel ' );
end

%% Compute the initial values

% Windows for computing initial values

if signal_len > (3*learn_period_len)
    
    learn_idx_1 = ceil( linspace( ...
        learn_period_len, ...
        signal_len - (2 * learn_period_len) + 1,  ...
        num_of_learning_segments ...
        ) );
    learn_idx_2 = learn_idx_1 + learn_period_len - 1;
    
else
    
    num_of_learning_segments = 1;
    learn_idx_1  = 1;
    learn_idx_2  = signal_len;
    
end

% initial "Xa"

mean_values = zeros( num_of_learning_segments, 1);

for s = 1 : num_of_learning_segments
    
    s_data = Z(learn_idx_1(s):learn_idx_2(s));
    s_sqi  = sqi(learn_idx_1(s):learn_idx_2(s));
    
    s_good_sqi = find(s_sqi >= sqi_threshold);
    
    if isempty(s_good_sqi)
        mean_values(s) = median( s_data );
    else
        mean_values(s) = median( s_data(s_good_sqi) );
    end

end

Xa(1) = median(mean_values);

Q(:)  = 1.0;
Pa(1) = 0.1;

% SQI corrected measurement noise variance 
% Li , Eq. 16
% R = R .* exp( ( 1 ./ ( sqi.^2) ) - 1 );

%% Estimation loop

for k = 1:signal_len
    
    %% Combining states
     Sigma = [];
     Innov_sum = 0;
     SQIcan_sum = 0;
     innov = 0 ;
     SQIcan = 0; 
     
    %calculating new X 
    for L = [1:length(Zs)]
        
        innov = Zs{L}(k) - Xa(k);
        
        SQIcan = sqis{L}(k);
        
        if innov <= 0.5 
            innov = 0.5;          
        end 
            
        Sigma(L) = (SQIcan/innov)^2 ; 
               
        %summing inovations and SQI 
        Innov_sum = Innov_sum + innov;
        
        SQIcan_sum = SQIcan_sum + SQIcan;
        
    end 
    
    % finding Xa components 
    Xa_sum = 0 ;
    R_sum = 0 ; 
    
    for L = [1:length(Zs)]
        
        XL = Zs{L}(k,2);
        
        summ = sum(Sigma);
        
        sigsig = Sigma(L);
        
        Xa_sum = Xa_sum + ( XL* (Sigma(L)/ sum(Sigma) )); %* (Innov_sum - Zs{L}(k) - Xa(k))/Innov_sum ;
        
        R_sum = R_sum + (sqis{L}(k,1) *  (Sigma(L)/ sum(Sigma) )); %* (Innov_sum - Zs{L}(k) - Xa(k))/Innov_sum ;
        
      
    end 

    %R_sum is the new SQI 
    R(k) = R(k) * exp( ( 1 ./ (R_sum^2) ) - 1 );
    
    Z(k) = Xa_sum ; 
    
    %% Measurement update (Correction)
    
    K(k) = Pa(k) / ( Pa(k) + R(k) );  % Bishop, Eq. 1.11
    
    r(k) = Z(k) - Xa(k);              % The innovation
    X(k) = Xa(k) + K(k) * r(k);       % Bishop, Eq. 1.12
    P(k) = ( 1 - K(k) ) * Pa(k);      % Bishop, Eq. 1.13
    
    % Time update (Prediction)
    
    Xa(k+1) = X(k);                   % Bishop, Eq. 1.9
    Pa(k+1) = P(k) + Q(k);            % Bishop, Eq. 1.10
    
end

% Correct array size, due to the "Time update" process for the last sample
Xa = Xa(1:signal_len);
Pa = Pa(1:signal_len);

%% Compute sqi by simply clipping the error covariance

% p_std = abs(P); 
% X_sqi = sqi;
% 
% X_sqi( isnan(X) ) = 0;
% X_sqi( isnan(X_sqi) ) = 0;
% 
% X_sqi( p_std > error_threshold) = 0;

%% Compute sqi by scaling: 2 * Aposteriori estimate error std

old_min_value    = min_error_threshold;
old_max_value    = max_error_threshold;
new_min_value    = 0;
new_max_value    = 1;
clip_values_flag = true;

P_2x_std = 2 * sqrt(abs(P));

p_scaled = pe.sigproc.change_number_range( ...
    P_2x_std, old_min_value, old_max_value, ...
    new_min_value, new_max_value, ...
    clip_values_flag, plot_detail_flag ...
    );

X_sqi = (p_scaled .* -1) + new_max_value;

%% Summary plot

if plot_summary_flag
    
    if isempty(t)
        t = (0:signal_len-1)' / fs;
    end
    
    p_rows = 6;
    p_cols = 1;
    p_row  = 0;
    ts_axh = [];
    
    figure;
    
    
    % Input data
    
    p_row = p_row + 1;
    ts_axh(p_row) = subplot(p_rows, p_cols, p_row );
    
    plot( t, Z);
    hold on;
    plot( t, X, 'r' );
    legend( 'Input', 'Kalman' );
    ylabel( 'Data' );
    
    
    % Input SQI
    
    p_row = p_row + 1;
    ts_axh(p_row) = subplot(p_rows, p_cols, p_row );
    
    plot( t, sqi);
    hold on;
    plot( t, X_sqi, 'r' );
    legend('Input', 'Output');
    ylabel( 'Input SQI' );
    set(gca,'YLim',[-0.1 1.1]);
    
    
    % Residual "r"
    
    p_row = p_row + 1;
    ts_axh(p_row) = subplot(p_rows, p_cols, p_row );
    
    plot( t, r );
    ylabel( {'Innovation' , 'r'} );
    
    % P
    
    p_row = p_row + 1;
    ts_axh(p_row) = subplot(p_rows, p_cols, p_row );
    
    plot( t, P );
    ylabel( 'P' );
    
    % Plot the predictions +/- 2 std of the error covariances
    
    p_row = p_row + 1;
    ts_axh(p_row) = subplot(p_rows, p_cols, p_row );
    
    up_limit  = X + P_2x_std;
    low_limit = X - P_2x_std;
    
    time_points  = [t;        flipdim(t,1)];
    value_points = [up_limit; flipdim(low_limit,1)];
    background_colour = [7 7 7]/8;
    
    min_value_points = min(value_points);
    max_value_points = max(value_points);
    
    fill(time_points, value_points, background_colour);
    hold on;
    plot(t, X , 'b', 'LineWidth', 2);
    ylabel( {'Error', 'covariance P'} );
    legend( '+/- 2*sqrt(P)', 'X');
    set(gca,'YLim',[min_value_points  max_value_points]);
    
    % Kalman good estimates
    
    p_row = p_row + 1;
    ts_axh(p_row) = subplot(p_rows, p_cols, p_row );
    
    x_good = X;
    x_good( X_sqi == 0) = NaN;
    
    min_y = nanmin(x_good)
    max_y = nanmax(x_good)
    
    plot( t, x_good );
    ylabel( { 'Kalman', 'good' } );
    set(gca,'YLim',[min_y  max_y]);
    
    xlabel( 'Time (seconds)'  );
    linkaxes(ts_axh,  'x');
    set( gca, 'XLim', [ t(1)  t(end) ] );
    
end

end

