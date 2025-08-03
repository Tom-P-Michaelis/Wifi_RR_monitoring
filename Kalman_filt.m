%%
% File    : Kalman_filt.m    
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
% DESCRIPTON:
% ----------
%
%  Kalman_filt 
%
%
% INPUT
% -----
%
% A: The time series predictions to be filtered. 
% 
% b: Vector which is matched to a.
% 
% FsRes: Frequency of the outputs.
%
% OUTPUT
% -----
% 
% 
% Lag: he time, in seconds, that vector b should be  
%  "advanced by" to optimally match the vectors


function Time_series = Kalman_filt(A,SQI)

%% Initializing
Kalman_A = [];

%initialising 
Q = 1;
R = 1;

%previous input as first input to start 
Xt1 = A(1);
Pt1 = 0; 

%% looping
for n = [1:length(A)]
    
    %current estimate as Zt 
    Zt = A(n);
     
    %prediction stage 
    X_t = Xt1; 
    
    P_t = Pt1 + Q; 

    %Update stage 
    
    Fac = 1/(SQI(n)^2) -1 ;
    
    R = 1*exp( Fac ) ;
    
    Kt = P_t*(P_t + R)^(-1);
    
    Xt = X_t + Kt*(Zt - X_t);
    
    Pt = (1 - Kt)*P_t;
    
    
    %% Storing     
    %posterior state estimates 
    
    Kalman_A(n) = Xt; 
    
    %storing for next prediction 
    Xt1 = Xt ;
    Pt1 = Pt;

end     

%assignment
Time_series = Kalman_A;

end 