%%
% File    : PCA_Filter.m    
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
%  PCA_Filter takes a windowed time series and seperates into principle
%  components. PCA_Filter then attemps to detect the compenent that is most
%  representative of the subjects breathing. This is done via computing a
%  score based on the SNR, the peak location and the presence of
%  Bimodality. The same process will be repeated for the HR.
%
%
% INPUT
% -----
%
% CSI: The windowed CSI time series on which PCA is performed. 
% The CSI contains 256 subcarriers to be reduced. 
% 
% Fs: The estimated sampling frequency of the the CSI. 
% 
% Time: The time of the windowed CSI - Used for plotting 
%
% OUTPUT
% -----
% 
% 
% RR: The principle component of the subject breathing as a time series
% vector. 
%
% HR: The principle component of the subject's cardiac cycle as a time series
% vector. 



function [Candidates, Index_Combined, Weights, Index_store]  = PCA_Candidates(CSI, Fs, Time)


     RR = [];
     RR_SNR = [];
     PCA_Peak =[] ;
     
     %% Detrending CSI 
     
     CSI_dt = detrend(CSI);
     CSI_dt = CSI_dt - mean(CSI_dt);
     CSI = CSI_dt;

    %% PCA decomposition 
    
   [~, scores, Weights] = pca(CSI);
    
    %% Hann Windowing for FFT 
    
    HW = hann(length(scores));    
    scores_hann = scores .* HW ; 

    %% Fourier transform of scores
    
    %fourier transform
    Length = length(scores);  
    
    Y_hann = fft(scores_hann);

    %remove sum term 
    FFT_scores = abs(Y_hann(1:floor(Length/2),:));
    Callibrate = (0:1:length(FFT_scores)-1)*Fs/(Length+1); 
    
    
    %% Analysing candidate PCs
    

    % Looping PCs and assigning a score
    
    % initializing 
    Index_store = [];  
    
    for n = [1:length(FFT_scores(1,:))/2] 
     
        %% Selecting data 
        
        % Selecting Data FFT
        FFT = FFT_scores(:,n);
        
        % Selectcting the variance of the component
        Var = Weights(n); 
        
         %% Mitigating 1/f noise 
%         
%         FFT_1F = FFT;
%         
%         Noise_Rep = ones(length(FFT_1F),1);
%         
        % Taking valid signal range as between 0.08 and 0.60  
        S1 = 0.08;
        S2 = 0.6;
%         
%         [~, SI1] = min(abs(Callibrate - S1));
%         
%         % increasing range to account for spilling 
%         
%         [~, SI2] = min(abs(Callibrate - S2));   
%         
%         Noise_Rep(SI1,1) = 1.5;
%         
%         Noise_Rep(SI1:SI2,1) = linspace(1.5, 1, length( Noise_Rep(SI1:SI2,1)));
%         
%         Noise_Rep(1:SI1,1) = linspace(1.5, 1.5, length(Noise_Rep(1:SI1,1)));
%                 
%        % combining 
%        
%        FFT_1F = FFT_1F ./ Noise_Rep;
%        
%        FFT = FFT_1F;

 %% Analysing VAR 
        
        VAR = Weights(n)/max(Weights); 
        
        % Removing very low Varience components 
        
        if  n <= 35      
            Var_index = 1 ;          
        else            
            Var_index = 0;  
        end 
        
        %% Analysing FFT peak 
        
        % ensuring peak is in the RR band 
        
        [~, Peak] = max(FFT);
        
        % using the previously defined criteria
        
        index_peak = 0;
        
        % if within breathing rate range
        if Callibrate(Peak) >= S1 && Callibrate(Peak) <= S2           
            index_peak = 1;             
        end 
        
        %% checkpointing 
        
        if index_peak == 1 && Var_index == 1 
                      
          
        %% Analysing SPI
        
        %creating time series 
        t = (0:length(scores(:,n))-1)' ./ Fs;
        
        % calculating spectral purity, a low SPI is a pure frequency domain
        [sp_index, ~, ~, ~] = pe.sigproc.spectral_purity(    ...
         scores(:,n)        , ...
        Fs                  , ...
        t                   , ...
        false                 ...
        );
    
        % inversing SPI so high SPI is a good signal 

%          sp_index = sp_index * (-1) + 1 ; 
%           
         sp_index = 1;
        %% Analysing SNR 
        
        % Taking SNR as the area in the valid range over total area. This
        % is equivalent to taking the average in the valid range
        
        Noise  = mean(FFT);
        
        % Taking valid signal range as between 0.08 and 0.60  
        S1 = 0.05;
        S2 = 0.65;
        
        [~, SI1] = min(abs(Callibrate - S1 - 0.04));
        
        % increasing range to account for spilling 
        
        [~, SI2] = min(abs(Callibrate - S2  +  0.2));      
        
        Signal = mean(FFT(SI1:SI2) );
        
        SNR = Signal/Noise ;
        
        
        else 
            
            SNR = 0 ;
            sp_index = 0; 
            
        end 
                      
        %% Storing Scores 
        
        Index_store(n,:) = [sp_index, SNR, Var_index, index_peak]; 
        
        %% Plotting scores and their PCs 
        
%         plot(Callibrate, FFT)
%         
%         pause
        
        
    end 
    
    %% Normalising Scores 
    
%     Index_store = [ normalize(Index_store(:,1),'range'), ... 
%                     normalize(Index_store(:,2),'range'), ... 
%                     normalize(Index_store(:,3),'range'), ... 
%                     normalize(Index_store(:,4),'range') ...
%                     ]; 

    %% Combining scores 
    
     Index_Combined = Index_store(:,1) .* Index_store(:,2) ...
         .* Index_store(:,3) .* Index_store(:,4);
    
    %% Candidate selection
    
    % finding x largest indexes 
    [Index , PCA_num] = maxk(Index_Combined, 5);
    
    % plotting best canditate
    
%     plot(Callibrate, FFT_scores(:,PCA_num(2)))
%     pause

    %% Filtering off Candidates 
    
    % initialising 
    Candidates = {};
    counter = 1 ;
    
    for m = [1:length(Index)]
        
        % arbitrary filter for i
        if Index(m) >= 0.01 
            
            % fisrt index refers to candidate rank, the second determines 
            Candidates{counter}{1} = [scores(:,PCA_num(m))];
            
            Candidates{counter}{2} = [Index(m)];
            
            Candidates{counter}{3} = [Index_store(PCA_num(m), :)];
       
            counter = counter + 1; 
            
      
        end 
        
    end 
    
    %% If empty array of candidates 
    
    if length(Candidates) <= 0
        
            Candidates{1}{1} = [scores(:,1)];
            Candidates{1}{2} = [0];
            Candidates{1}{3} = zeros(1,4);
    end 
    
    % Plotting Candidates 
     
%     plot([1:length(Candidates{1,1}{1,1})], Candidates{1,1}{1,1})
%     
%     size(Candidates{1,1}{1,1})
%     
%     
%     pause

%% Plotting 

% plotting first 10 PCs FFTs and their scores 
% 
% tiledlayout(11,1)
% 
% for n = [1:11] 
%     
%     nexttile 
%     plot(Callibrate, FFT_scores(:,n))
%     axis([0, Callibrate(length(Callibrate))/12, 0 , max(FFT_scores(:,n))*1.1])
%     
% end 


% 
% pause 



    
end 

