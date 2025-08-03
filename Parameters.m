%% Controlling Parameters

% Choosing which subjects to analyse. ["all", "b"] selects respiratory data
% from each participant.

 Subjects = ["all", "b"];
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

%%  Taking the HR 

for Participant = [1: length(Subjects)] 
    
    [CSI, Time_store, ~, ~] = Data_reader_V5(Subjects(Participant),"PCAP");
    
    
end 




