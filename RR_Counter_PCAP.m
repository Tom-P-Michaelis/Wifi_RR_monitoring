

function Predictions = RR_Counter_PCAP(Airflow, Fs, Time)
%%% function takes the CSV file that contains airflow data and the sampling
%%% frequency of the stowood for this file (currently 32Hz). It follows an
%%% algorithm proposed by Shaffer et al:
%%% (https://link.springer.com/content/pdf/10.1007/s10439-007-9428-1.pdf)
%%% and outputs the predicted respiritory rate as well as their
%%% corresponding time stamps.

%%% Still need to add - warning if failure 
%%% - how to handle saturates (Flat) airflow peaks
    
    Air_window_filt = Airflow;
    
    %% peak detection 
    
    
    Counter_maxs = 1 ;
    Maxs = [];
     for n = [2:length(Air_window_filt)-1]
         
         %finding if n is a stationary point 
         Turning = (Air_window_filt(n)-Air_window_filt(n-1))*(Air_window_filt(n+1)-Air_window_filt(n));
          
         %Storing stationary points
         if Turning < 0 
                       
             Maxs(Counter_maxs,1) = n ;
             Maxs(Counter_maxs,2) = Air_window_filt(n) ;
             Counter_maxs = Counter_maxs + 1 ;        
             
         end 
     end 
     
    %% Filtering off peaks & troughs 
    
    
    for l = [1:length(Maxs)-1]
        %difference of adjacent peaks 
        Maxs(l,3) = ((Maxs(l,2) - Maxs(l+1,2) )^2)^0.5;       
    end 
    
    %thresholding   
    Third_q = prctile(Maxs(:,3), 75);
    Min = 0;
    %Removing peak/trough if less than 0.3 of the given threshold
    while Min < Third_q*0.3 
        
       Min = min(Maxs(1:length(Maxs)-1,3)) ;
       
       P= 1;
       
       if  Min < Third_q*0.3 
       while P < length(Maxs)
           if Maxs(P,3) == Min           
               Maxs(P,:) = []; 
               Maxs(P,:) = []; 
               P=P-2;
           end 
           
           if P==-1
               P=0;
           end
           P=P+1;
       end 
       end 
       
       %recalculating differences
       for l = [1:length(Maxs)-1]
       %difference of adjacent peaks 
       Maxs(l,3) = (( Maxs(l,2) - Maxs(l+1,2) )^2)^0.5;  
       
    end 
    
    end 
     
    %% Extracting RR rate from peaks/troughs
    
    Start_time = Maxs(1,1)/Fs;
    
    %if even number of peaks dont take final peak
    if mod(length(Maxs(:,1)),2) == 0 
        Maxs(length(Maxs),:) = [] ;
    end 
    
    End_time = Maxs(length(Maxs),1)/Fs;
    
    %peaks in section
     Peak_no = (length(Maxs)-1)/2 ;
     
     %calculating RR 
     RR = 60*Peak_no/(End_time - Start_time);
    
    %% plotting 
    
%     if Time >= 902 && RR <= 20 
%         
%         
%    delete(gca)
%         
%     hold on 
%     
%     plot([1:length(Air_window_filt)]./Fs,Air_window_filt)
%     scatter(Maxs(:,1)./Fs, Maxs(:,2))

% xlabel("Time(s)")
% ylabel("CSI Value")
% 
%     
%     hold off
%     
%     Output = RR
%     
%     Time = Time
%     
%     pause
%     
%     
%     end 
    
%% storing RR and their time stamps 


Predictions = RR;
end 
