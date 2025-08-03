%%
% File    : CSI_predictor_V5.m    
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
%  CSI_predictor_V5 predicts the RR from CSI data. CSI_predictor_V5 calls
%  the function RR_Predict_WIFI_V5 which performs the prediction. 
%
%
% INPUT
% -----
%
% Subjects: The list of subject numbers to be included in the analysis 
% 
% 
% 
% OUTPUT
% -----
% 
% 
% RR_CSI: A cell vector containing the results of the RR predictions:
% 
% ┌───────
% │ RR_CSI    ├───┬──── ─── ───
% └───────      │
%                    │
%                    ▼
%                   P001b ─────────  Kalman ── ──
%                    │ │
%                    │ │
%                    │ └───────────── Median ── ──
%                    │
%                    │
%                    └─────────────── Unfiltered ──────►  Timestamp
% 
%                                               │    │
%                                               │    └────────────►  CSI Estimate1
%                                               │
%                                               │                                           
%                                               └─────────────────►  CSI Estimate2
% 
% _________________
%
% Improvement Aims from previous versions 
% -----
% 
% > Combine for all subjects  
% > In function filtering [DONE]
% > Clean & Comment
% > 

function [RR_CSI, SQI_out]  = CSI_predictor_V6(Subjects) 

%% Looping over all participants

% Initialising RR_CSI 

RR_CSI = {};
SQI_out = {};

for Participant = [1: length(Subjects)] 
    
    
    %% Reading Subject data 

    [CSI, Time_store, ~, ~] = Data_reader_V5(Subjects(Participant),"PCAP");

    CSI = abs(CSI); 

    %% Making Baseline RR Predictions

    [PCAP_predict, Linreg_predict, PCA_predict, SQI_store]  = RR_Predict_WIFI_V6(CSI, Time_store) ;
    
    %% Creating Results Table using best values
    

    
    SQI_Best = [] ;
    PCAP_predict_best = [];
    Linreg_predict_best = [];
    
    %for every prediction 
    for pred = [1:length(PCAP_predict{1})]
    
    SQI_temp = 0 ;
    SQI_temp_max = 0 ;
    % for each component
    for can = [1:length(PCAP_predict)]
        
        
        
        SQI_temp = SQI_store{can}(pred,:);
        
        if SQI_temp >= SQI_temp_max
            
            can_max = can;          
            
            SQI_temp_max = SQI_temp;
            
        end 
        
    end
    
    PCAP_predict_best(pred,:) = PCAP_predict{can_max}(pred,:);
    Linreg_predict_best(pred,:) = Linreg_predict{can_max}(pred,:);
    SQI_Best(pred) = SQI_store{can_max}(pred);
    
    end 
    
    %transposing
    SQI_Best = SQI_Best.';
    
    Results = [PCAP_predict_best(:,1), PCAP_predict_best(:,2), Linreg_predict_best(:,2), PCA_predict(:,2)];

    %% Median Filtering 

      %Median filter on  PCAP

      Results_median = [Results(:,1), medfilt1(Results(:,2) , 15,  'truncate') ...
            , medfilt1(Results(:,3) , 15,  'truncate') ,  medfilt1(Results(:,4) , 15,  'truncate')];

    %% Kalman Filtering 1D   
    
    % selecting best SQI data
    
    SQI = SQI_Best;
    
    % detrending the SQI
    SQI = detrend(SQI);
    
    % Normalising SQI 
    SQI = normalize(SQI,'range');
    
    T_fs = 1; 
    
    % applying kalman filters 
    plot_summary_flag = true;
    plot_detail_flag = true;
    
    min_error_threshold = 0.1;
    max_error_threshold = 3.3;

    
       [  CSI_Kalm , ...
            X_sqi  , ...
            r      , ...
            P        ...
         ] = ...
    pe.sigproc.kf.kalman_1d( ...
            Results(:,2)         , ...
            SQI                  , ...
            T_fs                  , ...
            min_error_threshold , ...
            max_error_threshold , ...
            plot_summary_flag   , ...
            plot_detail_flag    ...
        );
        
        [  Linreg_Kalm , ...
            ~  , ...
            r      , ...
            P        ...
         ] = ...
    pe.sigproc.kf.kalman_1d( ...
            Results(:,3)         , ...
            SQI                , ...
            T_fs                  , ...
            min_error_threshold , ...
            max_error_threshold , ...
            false   , ...
            false      );
        
    % saving results of kalman 
    Results_Kalman = [Results(:,1) , CSI_Kalm, Linreg_Kalm, Results(:,3)];
     
    %% Kalman Fusion bad 
    
%     % normalize SQI 
%     SQI_store_temp = [];
%     %convert to matrix
%     for N = [1:length(SQI_store)]          
%         SQI_store_temp(:,N) = SQI_store{N};   
%     end
%     % Normalize
%     
%     normA = SQI_store_temp - min(SQI_store_temp);
%     normA = normA ./ max(normA) ;
%     SQI_store_temp = normA;
%     
%     SQI_store_temp2 = {};
%     % return to original format 
%     for N = [1:length(SQI_store)]          
%         SQI_store_temp2{N} = SQI_store_temp(:,N);   
%     end
%     
%    SQI_store = SQI_store_temp2;
%     
%            [  CSI_Kalm , ...
%             X_sqi  , ...
%             r      , ...
%             P        ...
%          ] = ...
%     kalman_1d_fuse( ...
%             PCAP_predict         , ...
%             SQI_store                 , ...
%             T_fs                  , ...
%             min_error_threshold , ...
%             max_error_threshold , ...
%             plot_summary_flag   , ...
%             plot_detail_flag    ...
%         );
%         
%         [  Linreg_Kalm , ...
%             X_sqi  , ...
%             r      , ...
%             P        ...
%          ] = ...
%    kalman_1d_fuse( ...
%            Linreg_predict          , ...
%             SQI_store           , ...
%             T_fs                  , ...
%             min_error_threshold , ...
%             max_error_threshold , ...
%             false   , ...
%             false      );
%         
%     Results_Kalman = [Results(:,1) , CSI_Kalm, Linreg_Kalm, Results(:,3)];
    
    
    %% Kalman fusion good 
    
    % normalize SQI 
    SQI_store_temp = [];
    %convert to matrix
    for N = [1:length(SQI_store)]          
        SQI_store_temp(:,N) = SQI_store{N}; 
    end
    SQI_fusion = SQI_store_temp;
    
   SQI_fusion1 = SQI_fusion;
   % finding non 0 min SQI 
   minval = inf ;
   for i = [1:length(SQI_fusion1(:,1))]
       for j =  [1:length(SQI_fusion1(1,:))]
           
           val = SQI_fusion1(i,j);
           
           if val ~= 0 && val < minval
           minval = val;
           end           
       end 
   end 
    
    % Normalize step 1 
    for i = [1:length(SQI_fusion(:,1))]
       for j =  [1:length(SQI_fusion(1,:))]
           
           val = SQI_fusion1(i,j);
           
           if val ~= 0
           SQI_fusion1(i,j) = SQI_fusion1(i,j) - minval;
           end           
       end 
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % finding new maximum 
%    maxval = 0.001; 
%    for i = [1:length(SQI_fusion(:,1))]
%        for j =  [1:length(SQI_fusion(1,:))]
%            
%            val = SQI_fusion(i,j);
%            
%            if  val > maxval
%            maxval = val;
%            end           
%        end 
%    end 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

maxval = max(prctile(SQI_fusion1, 60), [], 'all');

if maxval == 0
    maxval = max(prctile(SQI_fusion1, 70), [], 'all');
end

if maxval == 0
    maxval = max(prctile(SQI_fusion1, 80), [], 'all');
end
   
    SQI_fusion1 = SQI_fusion1 ./ maxval;
   
    % ensuring nothing is above 0.9 
     for i = [1:length(SQI_fusion(:,1))]
       for j =  [1:length(SQI_fusion(1,:))]
           
           val = SQI_fusion1(i,j);
           
           if val >= 0.9
           SQI_fusion1(i,j) = 0.9;
           end           
       end 
    end 
    
    SQI_fusion = SQI_fusion1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     normA = SQI_store_temp - min(SQI_store_temp);
%     normA = normA ./ max(normA) ;
%     SQI_store_temp = normA;
%     
%     SQI_store_temp2 = {};
%     % return to original format 
%     for N = [1:length(SQI_store)]          
%         SQI_store_temp2{N} = SQI_store_temp(:,N);   
%     end
%     

% putting time series into matrix 

Time_series = [] ;
for n = [1:length(PCAP_predict)]
    
    Time_series(:,n) = PCAP_predict{n}(:,2);
    
end 


% Kalman Fusion 
    [ ...
            CSI_Kalm   , ...
            Z_out_sqi    ...
         ] = ...
     pe.sigproc.kf.kalman_fusion( ...
            Time_series                   , ...
            SQI_fusion                 , ...
            T_fs                  , ...
            min_error_threshold , ...
            max_error_threshold , ...
            plot_summary_flag   , ...
            plot_detail_flag     ...
        );
    
 % Kalman Fusion 
 % putting time series into matrix 
Time_series = [] ;
for n = [1:length(PCAP_predict)]
    
    Time_series(:,n) = Linreg_predict{n}(:,2);
    
end 


    [ ...
            Linreg_Kalm   , ...
            Z_out_sqi    ...
         ] = ...
     pe.sigproc.kf.kalman_fusion( ...
            Time_series                   , ...
            SQI_fusion                 , ...
            T_fs                  , ...
            min_error_threshold , ...
            max_error_threshold , ...
            plot_summary_flag   , ...
            plot_detail_flag     ...
        );   
    
    
    % saving results of kalman 
     Results_Kalman = [Results(:,1) , CSI_Kalm, Linreg_Kalm, Results(:,3)];
    
    
     %% Storing Results
     
     RR_CSI{Participant} = {Results , Results_median, Results_Kalman};
%      SQI_out{Participant} = [Z_out_sqi,X_sqi,SQI_Best];

%% SQI messing



end 


end 