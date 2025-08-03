%% Introduction
%
% File    : Wifi_Monitoring_V5.m    
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
% Wifi_Monitoring_V5 is the main script of the library. It calls all other
% functions in order to predict vital sign data from experimental results. 
%
% INPUT
% -----
%
% Wifi_Monitoring_V5 uses experimental data from Wi-Fi Vitals Monitoring 
% Dataset. This contains CSI and Stowood timeseries data. 
%
% OUTPUT
% -----
% Wifi_Monitoring_V5 outputs predictions in the form of the cell array
% Experimental_predictions. This can be saved to a CSV file. 
% 
% Wifi_Monitoring_V5 can output graphical and quantative analysis of the
% vital sign predictions. 
% ________________________________________________________________________
%
% Improvement Aims from V4 
% -----
% 
% > Implement ICA 
% > Extend to part a of experiment 
% > Calculate when data is not valid 
% > Fix SQI problems
% > Improve peak counting 
% > Re-order filtering and matching [DONE] 
% > Saving predictions optionally
% > Rearrange to avoid loops in the main document
%
% ________________________________________________________________________


%% Data Cleaning
% % 
 clc  
  clear all 

%% Path Maintenence 

addpath('C:\Users\tommi\OneDrive\Documents\4th Year Eng\4YP\Matlab\packages')

%% Controlling Parameters

% Choosing which subjects to analyse. ["all", "b"] selects respiratory data
% from each participant.

 Subjects = ["all", "b"];
%   Subjects = ["P002b"];
%  Subjects = ["all a and b"];

% choosing whether to save the predictions to a csv file
Save = "false";

%% Interpretting Subjects 

% If all subjects are chosen - return string array with every participants
% data

% dummy variable subs 
  Subs = [];    
  
if Subjects(1) == "all a and b" 
    
    for suffix = ["a", "b"] 
        
       for n = [1:9]
        
        Sub = strcat("P00" , num2str(n) , suffix);   
        Subs = [Subs ; Sub];
            
       end     
      for n = [10:15]

        
        Sub = strcat("P0" , num2str(n) , suffix);       
        Subs = [Subs ; Sub];
        
      end 
      Subjects = Subs;    
    end 
    
end   

if Subjects(1) == "all"   
  
    for n = [1:9]
        
        Sub = strcat("P00" , num2str(n) , Subjects(2));   
        Subs = [Subs ; Sub];
            
    end     
    for n = [10:15]

        
        Sub = strcat("P0" , num2str(n) , Subjects(2));       
        Subs = [Subs ; Sub];
            
    end    
    Subjects = Subs;    
end

%% Stowood Ground Truths 

Stowood_airflow_fs = 32 ;
Stowood_ppg_fs = 256; 

clear RR_Stowood

[RR_Stowood] = RR_Counter_Stowood_V10(Stowood_airflow_fs, Subjects);

%HR_Stowood = RR_Counter_V5(Stowood_ppg_fs, Subjects);

%% CSI based estimates 
% 
  [RR_CSI, SQI_out]  = CSI_predictor_V6(Subjects);

%% Postprocessing Data 
[Experimental_predictions] = DATA_handling_V5(RR_CSI, RR_Stowood, Subjects, Save);

%% Analysing Outputs

% save('Workspace')
 load('Workspace')

% closing all previous graphs 
    close all

% Calculating metrics and plotting results for all participants combined
[MAE , Interval, RMSE, TIME, MAELOW, TIME_TOTAL, RRAV] = Outputs_total_V5("all", Experimental_predictions, ...
        "Matched", "Kalman", "CSI", ["Bland comps"], Subjects, false)

% [MAE , Interval] = Outputs_total_V5("6", Experimental_predictions, ...
%     "Matched", "Kalman", "CSI", ["Results"])
%%



