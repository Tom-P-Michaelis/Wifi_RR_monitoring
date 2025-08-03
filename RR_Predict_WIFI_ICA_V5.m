%%
% File    : RR_Predict_WIFI_V5.m    
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
%  RR_Predict_WIFI_V5.m processed a complete time series CSI data with multiple
%  subcarriers. It performs windowing, PCA selection, filtering, peak
%  counting and FFT estimation. The function then calculates the predicted
%  BrPM from the processed signal.
%
%
%
% INPUT
% -----
%
% CSI: The entire time series CSI dataset with all carriers included. 
%
% Time_store: Time of each data point collected. 
% 
% OUTPUT
% -----
% 
% PCAP_predict: Standard CSI predictions that use a combination of peak
% counting and FFT prediction. 
%
% PCAP_predict: Predictions that only use the FFT method.  
%
% PCA_predict: Predictions that only use the FFT peak of the chosen principle
% component - only to be used for assesment of PCA performance.  
%
% SQI: Estimated signal quality metric.
% ________________________________________________________________________
%
% Improvement Aims from previous versions 
% -----
% 
% > Improve peak count to detect failure 
% > Combine for all subjects [DONE]
% > Change to power based filter selection 
% > Introduce results flagging 
% > Ouptput SQI


function [PCAP_predict, Linreg_predict, PCA_predict, SQI]  = RR_Predict_WIFI_ICA_V5(CSI, Time_store) 

%% CSI data capture rate 
%estimating the frequency of data capture

%selecting the last time stamp that is not 0
Time = 0;
n=-1;
while Time == 0 
   n = n + 1 ;
   Time = Time_store(length(Time_store)-n,1);  
end

%Computing the average capture rate in Hz 
Cap_rate = Time_store(length(Time_store)-n,1)/(length(Time_store)-n); 
Fs = 1/Cap_rate ;

%% Windowing 

%seconds for window size in seconds
Window_size = 30;

%defining windowing time resolution in seconds 
Interval = 1;

%creating grid (assuming linear spacing of CSI) 
Time_grid = [0:Cap_rate:(length(CSI)-1)*Cap_rate];

%Prediction storer
Predictions = [];
Linreg_predict = [];
PCA_predict_store = [];

%SQI storer
SQI_store = [];

%creating counters
counter = 1; 
z = 1; 

for Start_point = [1 : Interval*Fs : length(CSI) - Window_size*Fs]
 
    CSI_slice = CSI(round(Start_point):round(Start_point+Window_size*Fs),:);
    
    
    
    %% Subcarrier Selection (PCA) 
    
    %defining time of window for graphical purposes. 
    
    if mod(Window_size*Fs, 2) == 0
        Timestamp = Time_grid(round(Start_point + Window_size*Fs*0.5));
    else
        Timestamp = Time_grid(round(Start_point + (Window_size*Fs+1)*0.5));
    end
    
    [RR_CSI, RR_SNR, PCA_Peak] = ICA_Filter(CSI_slice, Fs, Timestamp);
    
    
    %% Filtering
    %Creating lowpass filter RR 
    
    %Checking if lowpass already exists exists
    A = exist('Filt_LP_RR', 'var');    
    if A == 1  
    else 
    Filt_LP_RR = designfilt('lowpassfir', ...
                         'PassbandFrequency', 40/60, 'StopbandFrequency', 80/60, ...
                         'SampleRate', Fs );
    end
    
    %Highpass filter RR
    % Filt_HP_RR = designfilt('highpassfir', ...
    %                      'PassbandFrequency', 6/60, 'StopbandFrequency', 3/60, ...
    %                      'SampleRate', 50);

    %RR_Filt = ButterworthRR_10thOrderm(Fs);

    % IIR_filt = filter(RR_Filt, RR_CSI);



    %double pass lp filter
    FIR_filt = filtfilt(Filt_LP_RR, RR_CSI);

    RR = FIR_filt;


    % plotting filter effect 

    % plot([1:length(RR_CSI)], RR_CSI , [1:length(IIR_filt)], IIR_filt, [1:length(RR)] , RR, 'LineWidth',0.1, 'LineWidth',2, 'LineWidth',2)
    % xlabel('CSI Section Packet number')
    % ylabel('CSI Value')
    % legend('Raw CSI section', 'IIR 10th order Butterworth Bandpass'  , 'FIR Lowpass')
    % 
    % pause

   %% Predictions using FFT  

   RR_FFT = fft(RR);
    
   %plotting 
   FFT_display = abs(RR_FFT(1:floor(length(RR)/2),:));
   Callibrate = (0:1:length(FFT_display)-1)*Fs/(length(RR)+1); 
   
   %plotting the filetered RRs and their FFTs
   
%    delete(gca) 
%    
%    tiledlayout(2,1)
%    
%    nexttile 
%    plot([1:length(RR)]./50 , RR)
%    
%    nexttile
%    plot(Callibrate, FFT_display)
%    
%    pause
%     
   
   %using linear regression on bins 
  
    Bin = 2;
    Omega_RR = Linreg(RR_FFT,Fs,Bin);
    

    %% Filtering predictions using peak counting 
    
    Peak_counter_predict = RR_Counter_PCAP_V5(RR, Fs, Time);
    
    % If predictions from peak counting very far from FFT peak - suggests
    % that peak counting has failed
    
%     Peak_counter_predict
%     
%     PCA_Peak*60
%     
%     (Peak_counter_predict - PCA_Peak*60)/(PCA_Peak*60)
    
%     if abs(Peak_counter_predict - PCA_Peak*60)/(PCA_Peak*60) >= 0.5
%         
%         Peak_counter_predict = PCA_Peak*60;
%         
%     end 
    
%     pause 
    
    %% Combining prediction sources
    
    % Taking the average of peak counting and linreg methods.
    Pred = (Omega_RR + Peak_counter_predict)/2;


    %% Plotting 
    
%     %plot the CSI slice
%     plot(Time_grid(round(Start_point):round(Start_point+Window_size*Fs)) , CSI_slice(:,10) )
%     
%     pause 
%     
%     delete(gca)
 
    %% Storing Predictions 
    
    %middle of window 
    if mod(Window_size*Fs, 2) == 0
        Timestamp = Time_grid(round(Start_point + Window_size*Fs*0.5));
    else
        Timestamp = Time_grid(round(Start_point + (Window_size*Fs+1)*0.5));
    end 
        
    Predictions(z,:) = [Timestamp, Pred] ;
    
    Linreg_predict(z,:) = [Timestamp, Omega_RR] ;
    
    PCA_predict_store(z,:) = [Timestamp, PCA_Peak*60] ;
    
    

    %% Calculating SQI    
    %SQI should vary between 0 (bad signal quality) and 1 (perfect signal
    %quality)
    
    %SQI should be lower if FFT based method and peak counting do not agree
    
    Agree_metric = abs((Omega_RR - Peak_counter_predict)) ;
   
    %not allowing inf values 
    if Agree_metric <= 0.5 && Agree_metric >= -0.5
        Agree_metric = 0.5;
    end 
    
    SQI_agree = Peak_counter_predict / Agree_metric ; 
    
    %SQI should be lower if outside of range that we expect 
    if Pred < 5 || Pred > 35 
        
        SQI_range = 0.1 ;
    else
        
        SQI_range = 1; 
    end
    
    %SQI should be larger if the PCA signal to noise ratio is higher 
    SQI_SNR = RR_SNR*Peak_counter_predict;
    
  
    
    %% Storing SQI 

    SQI_store(z,:) = [SQI_agree, SQI_range, SQI_SNR];
    
    z = z+1 ;
   
end 

%% outputting 

PCAP_predict = Predictions;
SQI = SQI_store;
PCA_predict = PCA_predict_store;

end 