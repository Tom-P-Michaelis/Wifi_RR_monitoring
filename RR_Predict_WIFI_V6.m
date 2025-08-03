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


function [PCAP_predict, Linreg_predict, PCA_predict, SQI]  = RR_Predict_WIFI_V6(CSI, Time_store) 

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

%% Clipping for speed 

%  CSI = CSI(100*Fs:200*Fs,:);

%% Windowing 

%seconds for window size in seconds
Window_size = 30;

%defining windowing time resolution in seconds 
Interval = 1;

%creating grid (assuming linear spacing of CSI) 
Time_grid = [0:Cap_rate:(length(CSI)-1)*Cap_rate];

%Prediction storer
Predictions = {};
Linreg_predict = {};
PCA_predict_store = [];

%SQI storer
SQI_store = {};

%creating counters
counter = 1; 
z = 1; 

Filter_config.method  = {'IIR'; 'IIR'};

for Start_point = [1 : Interval*Fs : length(CSI) - Window_size*Fs]
 
    CSI_slice = CSI(round(Start_point):round(Start_point+Window_size*Fs),:);
  
    
    %% Subcarrier Selection for Candidate Creation
    % Creating 20 candidate components for selection 
    
    
    %defining time of window for graphical purposes.     
    if mod(Window_size*Fs, 2) == 0
        Time = Time_grid(round(Start_point + Window_size*Fs*0.5));
    else
        Time = Time_grid(round(Start_point + (Window_size*Fs+1)*0.5));
    end
    
    % creating candidate components
    
    [Candidates, Index_Combined, Weights, Index_store] = PCA_Candidates(CSI_slice, Fs, Time); 
    
    %% Variance approach 
%     
VARs = [];
    for n = [1:length(CSI_slice(1,:))]
     subcarrier = CSI_slice(:,n);

     VARs(n) = var(subcarrier) ;
        
    end 

[val , num] = maxk(Index_Combined, 5);

for n = [1:length(Candidates)]
    
    
    Candidates{1,n}{1,1} = CSI_slice(:,num(n));

end 
    %% windowing
    
for PC = [1:length(Candidates)] 
        
        
    % selecting Candidate PC 
    RR_CSI = Candidates{1,PC}{1,1};
    
    
    % Selecting candidate SQI
%     SQI_Candidates = Candidates{1,PC}{1,2};
    
     %% SNR 
%    HW = hann(length(RR_CSI));    
%     RR_hann = RR_CSI .* HW ; 
% 
%     
%     %fourier transform
%     Length = length(RR_CSI);  
%     
%     Y_hann = abs(fft(RR_hann));
% 
%     %remove sum term 
%     FFT = abs(Y_hann(1:floor(Length/2),:));
%     Callibrate = (0:1:length(FFT)-1)*Fs/(Length+1); 
%     
%   
%       Noise  = mean(FFT);
%         
%         % Taking valid signal range as between 0.08 and 0.60  
%         S1 = 0.08;
%         S2 = 0.6;
%         
%         [~, SI1] = min(abs(Callibrate - S1 - 0.04));
%         
%         % increasing range to account for spilling 
%         
%         [~, SI2] = min(abs(Callibrate - S2  +  0.2));      
%         
%         Signal = mean(FFT(SI1:SI2) );
%         
%         SNR = Signal/Noise; 

% SNR = Index_store{PC}{
    
    %% Filtering IIR
    
    RR_CSI_N = RR_CSI; 
    
    % finding start value 
    Mem1 = median(RR_CSI(1:1*round(Fs)));
    
    % finding end value 
    Mem2 = median( RR_CSI(length(RR_CSI)-1*round(Fs):length(RR_CSI)) );
    
    % creating padding array that is 10 seconds long 
    
    Mem1 = Mem1 * ones(round(10*Fs),1)  ;
    Mem2 = Mem2 * ones(round(10*Fs),1)  ;
    
    % adding padding array 
    RR_CSI = [Mem1; RR_CSI; Mem2];
    
    [RR , Filter_config] = IIRFILT(RR_CSI, Fs, Filter_config, false, Time);
    
    % removing padding 
    
    RR = RR( round(10*Fs) : length(RR) - round(10*Fs) );
    
    RR_CSI = RR_CSI(round(10*Fs) : length(RR_CSI)- round(10*Fs));
    
    % No padding version 
    
    [RR_N , Filter_config] = IIRFILT(RR_CSI_N, Fs, Filter_config, false, 1);
    
    %% Plotting filter effect 
    
%     if Time >= 300
%         
%         t = tiledlayout(2,3)
%         
%         %plotting unfiltered
%         nexttile([1,2])
%         plot([1:length(RR_CSI)]./Fs, RR_CSI,  'LineWidth', 1.1)
%         
%         xlabel('Time (s)') 
%         ylabel('CSI') 
%         set(gca,'XMinorTick','on')
%         
%         axis([0, 30, min(RR_CSI)*1.1,  max(RR_CSI)*1.1])
%         
%         % plotting fft 
%         nexttile
%         
%         Y = fft(RR_CSI);
%         Y = abs(Y);
% 
%         Y_Half = (Y(1:floor(length(RR_CSI)/2),:));
%         Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(RR_CSI)+1); 
% 
%         plot(Callibrate, Y_Half,  'LineWidth', 1.1 )
%         
%         axis([0, 2 , 0 , max(Y_Half)*1.1])
%         
%         xlabel('Frequency (Hz)') 
%         ylabel('Power') 
%         set(gca,'XMinorTick','on')
%         
%         % plotting filtered
%         nexttile([1,2])
%         plot([1:length(RR_CSI)]./Fs, RR, 'r', 'LineWidth', 1.1)
%         
%         axis([0, 30, min(RR)*1.1,  max(RR)*1.1])
%         
%         xlabel('Time (s)') 
%         ylabel('CSI') 
%         set(gca,'XMinorTick','on')
%         
%         % next FFT 
%         
%         nexttile
%         
%         Y = fft(RR);
%         Y = abs(Y);
% 
%         Y_Half = (Y(1:floor(length(RR)/2),:));
%         Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(RR)+1); 
% 
%         plot(Callibrate, Y_Half, 'r',  'LineWidth', 1.1 )
%         
%         axis([0, 2 , 0 , max(Y_Half)*1.1])
%         
%         xlabel('Frequency (Hz)') 
%         ylabel('Power') 
%         set(gca,'XMinorTick','on')
% 
%         pause   
%         delete(gca)
%         
%        %% plotting Hann effect on filtered 
%        
%        
%                
%         t = tiledlayout(2,3)
%         
%         %plotting unfiltered
%         nexttile([1,2])
%         plot([1:length(RR)]./Fs, RR,  'LineWidth', 1.1)
%         
%         xlabel('Time (s)') 
%         ylabel('CSI') 
%         set(gca,'XMinorTick','on')
%         
%         axis([0, 30, min(RR)*1.1,  max(RR)*1.1])
%         
%         % plotting fft 
%         nexttile
%         
%         Y = fft(RR);
%         Y = abs(Y);
% 
%         Y_Half = (Y(1:floor(length(RR)/2),:));
%         Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(RR)+1); 
% 
%         plot(Callibrate, Y_Half,  'LineWidth', 1.1 )
%         
%         axis([0, 1 , 0 , max(Y_Half)*1.1])
%         
%         xlabel('Frequency (Hz)') 
%         ylabel('Power') 
%         set(gca,'XMinorTick','on')
%         
%         % plotting filtered
%         
%         RR = RR .* hann(length(RR));
% 
%         nexttile([1,2])
%         plot([1:length(RR_CSI)]./Fs, RR, 'r', 'LineWidth', 1.1)
%         
%         axis([0, 30, min(RR)*1.1,  max(RR)*1.1])
%         
%         xlabel('Time (s)') 
%         ylabel('CSI') 
%         set(gca,'XMinorTick','on')
%         
%         % next FFT 
%         
%         nexttile
%         
%         Y = fft(RR);
%         Y = abs(Y);
% 
%         Y_Half = (Y(1:floor(length(RR)/2),:));
%         Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(RR)+1); 
% 
%         plot(Callibrate, Y_Half, 'r',  'LineWidth', 1.1 )
%         
%         axis([0, 1 , 0 , max(Y_Half)*1.1])
%         
%         xlabel('Frequency (Hz)') 
%         ylabel('Power') 
%         set(gca,'XMinorTick','on')
%         
%         hold on 
%         
%         [peak, index] = max(Y_Half);
%         scatter(Callibrate(index),peak, 'k', 'LineWidth', 2) 
%         
%         legend('FFT Peak 0.2Hz (12 brpm)')
% 
%         pause   
%         delete(gca)
%        
% 
%     end 

   %% Predictions using FFT  

   RR_Hann = RR .* hann(length(RR));
   RR_FFT = fft(RR_Hann);
    
   %plotting 
   FFT_display = abs(RR_FFT(1:floor(length(RR)/2),:));
   Callibrate = (0:1:length(FFT_display)-1)*Fs/(length(RR)+1); 
   
   %plotting the filtered RRs and their FFTs
   
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


    
    %% taking FFT Peak as Omega 
    
    Y = RR_FFT ;
    len = length(Y);

    FFT_half= (Y(1:floor(len/2)));
    Callibrate = (0:1:length(FFT_half)-1)*Fs/(len+1); 


    FFT_RR = FFT_half;

    [~, freq] = max(FFT_half);

    %frequency is target frequency 
    Omega_RR = Callibrate(freq)*60;
    
      
%    Bin = 2;
%     Omega_RR = Linreg(RR_FFT,Fs,Bin);

    %% Filtering predictions using peak counting 
    
    Peak_counter_predict = RR_Counter_PCAP_V6(RR, Fs, Time);
    
    
    %% Combining prediction sources
    
    % Taking the average of peak counting and linreg methods.
    Pred = (Peak_counter_predict+Omega_RR)/2;


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
        
    Predictions{PC}(z,:) = [Timestamp, Pred] ;
    
    Linreg_predict{PC}(z,:) = [Timestamp, Omega_RR] ;
    
    % PCA peak should be nulled 
%     PCA_predict_store(z,:) = [Timestamp, PCA_Peak*60] ;
    PCA_predict_store(z,:) = [Timestamp, 0.5*60] ;
  
    %% Creating another SQI variable based on peak count vs linreg
        
    %SQI should be lower if FFT based method and peak counting do not agree
    
    Agree_metric = abs((Omega_RR - Peak_counter_predict)) ;
   
%     %not allowing inf values 
%     if Agree_metric <= 1 
%         Agree_metric = 1;
%     end 
%     
%     SQI_agree = 1 / Agree_metric ;

% Binarising SQI agree 

if Agree_metric >= 3
    SQI_agree = 0;
else 
    SQI_agree = 1;
end 
    
    %% Getting SPI SQI 
    
    SQI_Candidates = Candidates{PC}{1,3}(1,1);
    
    %% Recalculating SPI 
           
        %creating time series 
        t = (0:length(RR)-1)' ./ Fs;
        
        % calculating spectral purity, a low SPI is a pure frequency domain
        [sp_index, ~, ~, ~] = pe.sigproc.spectral_purity(    ...
         RR        , ...
        Fs                  , ...
        t                   , ...
        false                 ...
        );
    
    if Pred <= 11
        sp_index = sp_index*1.33;
    end 
        % inversing SPI so high SPI is a good signal 

%         sp_index = sp_index * (-1) + 1 ; 


%  sp_index = Candidates{PC}{1,3}(1,1);

 SNR = Candidates{PC}{1,3}(1,2);
     
    %% SQI Admin    
    
     SQI_store{PC}(z,:) = [SNR*SQI_agree*sp_index];
% SQI_store{PC}(z,:) = [SQI_Candidates*SQI_agree];


    end 
    z = z+1 ;
   
   
end 

%% Cell array admin 

normal_length = length(Predictions{1});

for n = [1:length(Predictions)]
    
    if length(Predictions{n})< normal_length
        
        Predictions{n} = [Predictions{n} ;  zeros(normal_length- length(Predictions{n}), 2)];
        Linreg_predict{n} = [Linreg_predict{n} ;  zeros(normal_length- length(Linreg_predict{n}), 2)];
        SQI_store{n} = [SQI_store{n} ;  zeros(normal_length- length(SQI_store{n}), 1)];
        
    end 
end 



%% outputting 

PCAP_predict = Predictions;
SQI = SQI_store;
PCA_predict = PCA_predict_store;

end 