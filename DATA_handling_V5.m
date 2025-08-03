% 
% 
% 
% ┌────────────┐
% │Experimental       │
% └─┬──────┬───┘
%    │          │
%    │          │
%    │          └──►
%    │
%    ▼
%    RR ───►
%    │
%    │
%    │
%    │
%    ▼
%   P001b
%   │   │
%   │   └───────────┐
%   ▼                       ▼
%   Matched              Unmatched ── ──
%   │ │  │
%   │ │  │
%   │ │  │
%   │ │  └─────── Kalman ── ──
%   │ │
%   │ │
%   │ └────────── Median ── ──
%   │
%   │
%   │
%   ▼
% Unfiltered ──────►  Timestamp
% │ │ │ │ │
% │ │ │ │ │
% │ │ │ │ └─────────►  Stowood GT
% │ │ │ │
% │ │ │ └───────────►  CSI predictions
% │ │ │
% │ │ └─────────────►  Linreg
% │ │
% │ └───────────────►  PCA predictions
% │
% └──────────────────► Metronome Target





function [Experimental_predictions] = DATA_handling_V5(RR_CSI, RR_Stowood, Subjects, Save)

%% Data Initialisation 
Experimental_predictions = {};

FsRes = 1 ;


%% Metronome Data Creation 

    %recreating the target RR with arbitrary start time 
    Metronome = [6 : 3 : 33];

    %timestamps with frequency FsRes 
    Metronome_Time = [1 : 1/FsRes: length(Metronome)*120];

    Metronome_target = zeros(1200,2);
    
    %initializing counter 
    counter = 1; 
    for n = Metronome_Time

        Metronome_target(n,:) = [Metronome_Time(n) , Metronome(counter)];

        if mod(n,120) == 0 
        counter = counter + 1;
        end 

    end 

    %Renaming 
    Metronome_predict = Metronome_target; 


   %% Looping over Participants
for Participant = [1: length(Subjects)] 
    
    %initialising lags
    Lag_PCAP = [] ;
    Stowood_predict = [];
    Linreg_predict = [];
    PCA_predict = [];
    
   %% Looping over Method
   for Method = [1:3]  
       
    
       %% DATA Extraction 
       
       Times = RR_CSI{Participant}{Method}(:,1);
       
       PCAP_predict = [Times, RR_CSI{Participant}{Method}(:,2)];
       
%        size(PCAP_predict)
       
       Linreg_predict = [Times, RR_CSI{Participant}{Method}(:,3)];
       
%        size(Linreg_predict)
       
       PCA_predict = [Times, RR_CSI{Participant}{Method}(:,4)];
       
%        size(PCA_predict)
       
       Stowood_predict = RR_Stowood{Participant}(:,:);
       
%        size(Stowood_predict)
       
       Metronome_prediction = Metronome_predict;
       
%        size(Metronome_prediction)
       
    %% Cutting endpoints

    %cutting the first period off of all results 
    Cut1 = 45 ;
    
    %cutting the last period off of all results 
    Cut2 = 30;

    % Removing data start 
    PCAP_predict([1:round(Cut1*FsRes)], :) = [];
    Stowood_predict([1:round(Cut1*FsRes)], :) = [];
    Metronome_prediction([1:round(Cut1*FsRes)], :) = [];
    Linreg_predict([1:round(Cut1*FsRes)], :) = [];
    PCA_predict([1:round(Cut1*FsRes)], :) = [];
    
%     size(PCA_predict)
    
    %cutting the end of the dataset 
    
    %cutting metronome first, 
    Metronome_prediction(length(Metronome_prediction)-[0:round(Cut2*FsRes)-1], :) = [];

    % cutting PCAP and Stowood to the same length
%    min(length(Metronome_prediction),length(PCAP_predict))
    
    PCAP_predict(length(Metronome_prediction) +1 : length(PCAP_predict) , :) = [];
    Stowood_predict(length(Metronome_prediction) +1  : length(Stowood_predict), : ) = [];
    Linreg_predict(length(Metronome_prediction) +1 : length(Linreg_predict), : ) = [];
    PCA_predict(length(Metronome_prediction) +1  :  length(PCA_predict), : ) = [];
 
    
    
    %% Matching 
    %Lags are the length of time, in seconds, that the results should be lagged
    %to match the stowood data. 

    %matching stowood and pcap 
    Lag_PCAP(Method) = Matching(Stowood_predict, PCAP_predict, FsRes);
  
    %matching stowood and Metronome  
    Lag_Metronome(Method) = Matching(Stowood_predict, Metronome_prediction, FsRes);
  
    %matching stowood and Linreg  
    Lag_Linreg(Method) = Matching(Stowood_predict, Linreg_predict, FsRes);
    
    %% Creating unmatched results by combining results so far 
 
    % array concatenation consistency.
        if length(PCA_predict)== length(Stowood_predict)
            
            Results_Unmatched = [Stowood_predict(:,1), Stowood_predict(:,2)...
                , PCAP_predict(:,2), Linreg_predict(:,2),...
                PCA_predict(:,2), Metronome_prediction([1:length(PCAP_predict)],2)];
            
        elseif length(PCA_predict)>= length(Stowood_predict)
            
            Results_Unmatched = [Stowood_predict(:,1), Stowood_predict(:,2)...
                , PCAP_predict([1:length(Stowood_predict)],2),...
                Linreg_predict([1:length(Stowood_predict)],2),...
                PCA_predict([1:length(Stowood_predict)],2),...
                Metronome_prediction([1:length(Stowood_predict)],2)];
        end 


    %% Combining Matched Results

    %apllying the lags and combining into Results table
    Results = [];
    N=1;
    
    % For all other filters use their lag   
    if Method ~= 5
    
    for n = [1:length(Stowood_predict)]

        if n + Lag_PCAP(Method) > 0 && n + Lag_PCAP(Method) <= length(PCAP_predict) ...
                && n + Lag_Metronome(Method) > 0 && ...
                n +  Lag_Metronome(Method) <= length(Metronome_prediction) && ...
                n + Lag_Linreg(Method) > 0 && n + Lag_Linreg(Method) <= length(Linreg_predict) ...
                &&  n + Lag_Linreg(Method) > 0 && n + Lag_Linreg(Method) <= length(PCA_predict) ...
                
            
        %N is the new timestamp from 0     
        Results(N,:) = [N , Stowood_predict(n,2), PCAP_predict(n+Lag_PCAP(Method),2) ,...
             Linreg_predict(n + Lag_Linreg(Method), 2), PCA_predict(n + Lag_Linreg(Method), 2), ...
             Metronome_prediction(n+ Lag_Metronome(Method),2)];
        
        N = N+1;
        end 

    end
    
    else
        display([Lag_PCAP(2), Lag_PCAP(3)])
        
     for n = [1:length(Stowood_predict)]

        if n + Lag_PCAP(2) > 0 && n + Lag_PCAP(2) <= length(PCAP_predict) ...
                && n + Lag_Metronome(2) > 0 && ...
                n +  Lag_Metronome(2) <= length(Metronome_prediction) && ...
                n + Lag_Linreg(2) > 0 && n + Lag_Linreg(2) <= length(Linreg_predict) ...
                &&  n + Lag_Linreg(2) > 0 && n + Lag_Linreg(2) <= length(PCA_predict) ...
                
            
        %N is the new timestamp from 0     
        Results(N,:) = [N , Stowood_predict(n,2), PCAP_predict(n+Lag_PCAP(2),2) ,...
             Linreg_predict(n + Lag_Linreg(2), 2), PCA_predict(n + Lag_Linreg(2), 2), ...
             Metronome_prediction(n+ Lag_Metronome(2),2)];
        
        N = N+1;
        end 


    end
        
    end 
    
    % Results is now a matrix containing timestamps of result (at 1HZ) in
    % column 1, Stowood RR predictions in column 2, peak counting results
    % in column 3, metronome target results in column 4 and linear
    % regression results in coloumn 5. 
    
    %% Flagging Unreliable Stowood Data
    
    %
    
    
    %% Storing Matched & Unmatched DATA 
    
    Experimental_predictions{Participant}{Method}{1} = Results;
    Experimental_predictions{Participant}{Method}{2} = Results_Unmatched;
    
    %% Saving Data 
    
    % if save is true 
    if Save == "True"
        save('Experimental_predictions.mat', 'Experimental_predictions');
    end 
    
   end 

end 
 
end 