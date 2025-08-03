%%
% File    : Outputs_total_V5.m    
% Author  : Tom Michaelis (tom.michaelis.tm@googlemail.com)
% Created : 12/04/2021
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
%  Outputs_total_V5.m calculates accuracy metrics between 2 different sources
%  over all participants.
%
%
% INPUT
% -----
%
% Participant: A string which represents the number for the participant
% data being used. "all" uses every available participant. "8" only uses
% the 8th participant. 
% 
% Experimental_predictions: The results in a cell array.
% 
% Matched: Dictates wether matched or unmatched data should be used.
% Options are "Matched" or "Unmatched".
%
% Filter: Dictates which filtering option to take data from.
% Options are "Unfiltered", "Median" or "Kalman".
%
% Source: Dictates which prediction source to take data from. Options are  
% "CSI" (standard), "Linreg", "Metronome" or "PCA".
%
% Plot: Decides which output plots to make. Can be a string vector
% containing a combination of: "Matching Results", "Results", "Bland",
% "HOR", "HOD" or "Interval HOR"
%
% OUTPUT
% -----
% 
% MAE: Mean absolute error for the chosen configuration 
%
% Interval: The percentage of results that are within 2BrPM of the target. 
%


function [MAE , Interval, RMSE, TIME, MAELOW, TIME_TOTAL, RRAV] = Outputs_total_V5(Participant, Experimental_predictions,Matched, Filter, Source, Plot, Subjects, movement) 

    
%% Calculating Metrics

% creating metrics
AE_PCAP = 0 ;
SE_PCAP = 0 ;
AE_LOW_PCAP = 0;

Noflag = 0 ;
Orig_Time = 0 ;

%initializing 
Low_Count = 0 ;
Good = 0;
m = 1 ;
mm = 1;
mm2 = 1;
mm3 = 1;
Differences = [];
Means = [];
Bland_Store = [];
Bland_Store1 = [];
Bland_Store2 = [];
Bland_Store3 = []; 
Grounds = [];
Ground = 2;
Index_mode = 0;


%% Choosing the array to take predictions from.

% Making sense of Matched, Filter and Source inputs

% Interpretting Matched
if Matched ==  "Matched" 
    Matched_index = 1;
elseif Matched ==  "Unmatched" 
    Matched_index = 2;
end 

% Interpretting Filter
if Filter == "Unfiltered"  
    Filter_index = 1;
elseif Filter == "Median"
    Filter_index = 2;    
elseif Filter == "Kalman" 
    Filter_index = 3;
end 

% Interpretting Source
if Source == "CSI"   
    Source_index = 3;
elseif Source == "Linreg"   
    Source_index = 4 ; 
elseif Source == "PCA"   
    Source_index = 5 ;    
elseif Source == "Metronome"   
    Source_index = 6 ;  
end 

%interpretting participant choice
if Participant == "all"
    cell = [1:length(Experimental_predictions)];
else 
    cell = str2num(Participant);   
end 
    
%% Looping over each participant 

% performing type b results analysis
    
for cell = cell
    
    %% interpretting subjects 

    A = Subjects(cell);

    B = convertStringsToChars(A);

    Type = B(length(B));
    
    
    if Type == 'b'
    %% 
	 
    %defining ground and estimates vectors
    Ground_truth_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,2);
    Estimate_value_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,Source_index);
    Time_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,1);
    Metronome_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,6);
    
    % initialising stores 
    Flag = zeros(length(Ground_truth_vector),1);
    Flag_Bloat = zeros(length(Ground_truth_vector),1);

    
    %% Looping over ground and estimate arrays to find bad Stowood predictions 
    
    for n = [1:length(Ground_truth_vector)]
        
    %Defining estimate, ground truth, time, and Metronome
    Ground_truth = Ground_truth_vector(n);
    Estimate_value = Estimate_value_vector(n);
    Metronome_value = Metronome_vector(n);
    
    %% Flagging Bad Stowood (No flagging if source is metronome) 
    
    if abs(Metronome_value - Ground_truth) >= 4 && Source ~= "Metronome"   
        
        % setting flag to true 
        Flag(n) = 1 ;
               
    end   
    
    end
    
    %% Bloating Flags so that adjacent values are also binned.
    Flags_Bloat = [];
    
    for r = [1:length(Flag)] 
        % adjacent value selection 
        for bin = [-5:5]
            % ensuring within index 
            if r+bin > 0 &&  r+bin <= length(Flag)
                %if adjacent flag is on
                if Flag(r+bin)
                    Flag_Bloat(r) = 1 ;
                end     
            end  
        end
    end 

    %Renaming 
    Flag = Flag_Bloat ;
    
    
for n = [1:length(Ground_truth_vector)]
    
    %Defining estimate, ground truth, time, and Metronome
    Ground_truth = Ground_truth_vector(n);
    Estimate_value = Estimate_value_vector(n);
    Metronome_value = Metronome_vector(n);
    
    %% Not Calculating values if metronome has failed
    Orig_Time = Orig_Time + 1 ;
    
    if Flag(n) == 1
              
    else 
    Noflag = Noflag + 1;
    %% Calculating Metrics 
    %Summing absolute error 
     
    AE_PCAP = AE_PCAP + abs(Ground_truth - Estimate_value);
    SE_PCAP = SE_PCAP + (Ground_truth - Estimate_value)^2;
    
    if Ground_truth < 12 
        AE_LOW_PCAP = AE_LOW_PCAP + abs(Ground_truth - Estimate_value);
        Low_Count = Low_Count + 1; 
    end 
    
    %counting "good" results as within 2
    if abs(Ground_truth -  Estimate_value) <= 2        
        Good = Good + 1;         
    end 
    
    Differences(m) = [Ground_truth -  Estimate_value];
    
    Means(m) = (Ground_truth + Estimate_value)/2;
    
    Bland_Store(m,:) = [Ground_truth, Estimate_value];
    
    Grounds(m) = Ground_truth;
    
    if Ground_truth <= 12 && Ground_truth > 2 && Estimate_value > 2
        Bland_Store1(mm,:) = [Ground_truth, Estimate_value];
        mm = mm + 1; 
    end 
    
    if Ground_truth > 12 && Ground_truth < 21 && Estimate_value > 2
        Bland_Store2(mm2,:) = [Ground_truth, Estimate_value];
        mm2 = mm2 + 1; 
    end 
    
    if Ground_truth >= 21 && Ground_truth <=35 && Estimate_value > 2
        Bland_Store3(mm3,:) = [Ground_truth, Estimate_value];
        mm3 = mm3 + 1; 
    end 
    
    
    % m defines how many iterations have been done
    m = m + 1;
    
    end 
    
end

    %% For HR experiment 
    
    elseif Type == 'a'
        
	 
    %defining ground and estimates vectors
    Ground_truth_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,2);
    Estimate_value_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,Source_index);
    Time_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,1);
    Metronome_vector = Experimental_predictions{cell}{Filter_index}{Matched_index}(:,6);
    
   
    
    
for n = [1:length(Ground_truth_vector)]
    
    %Defining estimate, ground truth, time, and Metronome
    Ground_truth = Ground_truth_vector(n);
    Estimate_value = Estimate_value_vector(n);
    Metronome_value = Metronome_vector(n);
    
    %% Not Calculating values if in bad timings
    
    % I already chop off 45 at the start and 30 at the end. 
    
    % from 1m to 3.75m, then 14m to 17m
    if movement == false && n >= 15 && n <= 225  || ...
        movement == false && n >= 840 &&  n <= 1020
        
    %% Calculating Metrics 
    %Summing absolute error 
    AE_PCAP = AE_PCAP + abs(Ground_truth - Estimate_value);
    SE_PCAP = SE_PCAP + (Ground_truth - Estimate_value)^2;
    
    %counting "good" results as within 2
    if abs(Ground_truth -  Estimate_value) <= 2        
        Good = Good + 1;         
    end 
    
    Differences(m) = [Ground_truth -  Estimate_value];
    
    Means(m) = (Ground_truth + Estimate_value)/2;
    
    Bland_Store(m,:) = [Ground_truth, Estimate_value];
    
    Grounds(m) = Ground_truth;
    
    % m defines how many iterations have been done
    m = m + 1;
     
%     Repeating if movement is on 

    elseif movement == true && n >= 360 && n <= 540  
        
    %% Calculating Metrics 
    %Summing absolute error 
    AE_PCAP = AE_PCAP + abs(Ground_truth - Estimate_value);
    SE_PCAP = SE_PCAP + (Ground_truth - Estimate_value)^2;
    
    %counting "good" results as within 2
    if abs(Ground_truth -  Estimate_value) <= 2        
        Good = Good + 1;         
    end 
    
    Differences(m) = [Ground_truth -  Estimate_value];
    
    Means(m) = (Ground_truth + Estimate_value)/2;
    
    Bland_Store(m,:) = [Ground_truth, Estimate_value];
    
    Grounds(m) = Ground_truth;
    
    
    % m defines how many iterations have been done
    m = m + 1;
     
    
   
    end
        
    end 
    end 
end
    
    

%% Calculating total metrics 

% percentage of results within 2 BrPM
Interval = (Good/(m-1))*100;
% Mean absolute error
MAE = AE_PCAP/(m-1);
RMSE = (SE_PCAP/(m-1))^0.5;
MAELOW = AE_LOW_PCAP/Low_Count ;

TIME = Noflag/60;
TIME_TOTAL = Orig_Time/60;

RRAV = std(Grounds);

%% Plotting

% plotting Bland Altman, correlation plot, histogram of differences and
% mean results

if ismember("Bland", Plot)  == 1
    
    %Bland Altman Plot 
    BlandAltman(Bland_Store(:,1) , Bland_Store(:,2) ,{'Reference Estimates (brpm)','CSI Estimates (brpm)'},'', [], 'baInfo', {'IQR'}, 'corrInfo', {'r2'}, 'showFitCI', 'on', 'diffValueMode', 'Absolute', "data1Mode", 'Truth',  'baStatsMode', 'Non-parametric' );
axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 15  % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(b)',  'fontweight', 'bold')

x_txt = min(axlim) - 85  % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(a)',  'fontweight', 'bold')
end 

if ismember("Bland low", Plot)  == 1
    
    %Bland Altman Plot 
    BlandAltman(Bland_Store1(:,1) , Bland_Store1(:,2) ,{'Reference Estimates (brpm)','CSI Estimates (brpm)'},'', [], 'baInfo', 'IQR', 'corrInfo', {'r2'}, 'showFitCI', 'on', 'diffValueMode', 'Absolute', "data1Mode", 'Truth',  'baStatsMode', 'Non-parametric' , 'axesLimits', [3,18], 'baYLimMode', 'auto');
    
axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 3.5  % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(b)',  'fontweight', 'bold')

x_txt = min(axlim) - 22  % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(a)',  'fontweight', 'bold')

end 

if ismember("Bland comps", Plot)  == 1
    
    %Bland Altman Plot 
    
    BlandAltman(Bland_Store1(:,1) , Bland_Store1(:,2) ,{'Reference Estimates (brpm)','CSI Estimates (brpm)'},'', [], 'baInfo', 'IQR', 'corrInfo', {'r2'}, 'showFitCI', 'on', 'diffValueMode', 'Absolute', "data1Mode", 'Truth',  'baStatsMode', 'Non-parametric' , 'axesLimits', [3,18], 'baYLimMode', 'Square');
    BlandAltman(Bland_Store2(:,1) , Bland_Store2(:,2) ,{'Reference Estimates (brpm)','CSI Estimates (brpm)'},'', [], 'baInfo', 'IQR', 'corrInfo', {'r2'}, 'showFitCI', 'on', 'diffValueMode', 'Absolute', "data1Mode", 'Truth',  'baStatsMode', 'Non-parametric' , 'axesLimits', [3,18], 'baYLimMode', 'Square');
    BlandAltman(Bland_Store3(:,1) , Bland_Store3(:,2) ,{'Reference Estimates (brpm)','CSI Estimates (brpm)'},'', [], 'baInfo', 'IQR', 'corrInfo', {'r2'}, 'showFitCI', 'on', 'diffValueMode', 'Absolute', "data1Mode", 'Truth',  'baStatsMode', 'Non-parametric' , 'axesLimits', [3,18], 'baYLimMode', 'Square');


end 

%Histogram plots 
if ismember("HOD", Plot)  == 1
    
    %Histogram of differences 
    figure
    histogram(Differences)
    title('Histogram of differences')
    xlabel('Reference RR - CSI RR (brpm)')
    ylabel('Estimates')
    axis([-12, 12, 0, 2700])
    
end 

%Histogram plots 
if ismember("HOR", Plot)  == 1
    
    
    %Histogram of results
    figure 
    H = histogram(Means, 10)
    xlabel('Mean of Reference RR and CSI RR')
    ylabel('Count')
    
    
    
    
    %Histogram of Stowood results
    
    figure 
  
    
    histogram(Grounds, 10)
    xlabel('Respiratory Rate (brpm)')
    ylabel('Estimates')
      
    set(gca,'XMinorTick','on')
    
    hold on     
    
    J = boxplot(Grounds, 'orientation', 'horizontal', 'Positions', (max(H.Values)*1.2),'Widths',max(H.Values)*0.1, 'Whisker',1 );
    
    set(J,{'linew'},{1})
    yticks([])
    
    % resetting y ticks 
    tik = [0:1000:max(H.Values)*1.2];
    yticks(tik)
    
    labs = {};
    for n = [1:length(tik)]
        
        labs{n} = char(num2str(tik(n)));      
    end 
    
    yticklabels(labs)
    
    axis([0,36,0, max(H.Values)*1.4]) 
    
    hold off 
    
end

% Histogram of Results with confidence intervals
if ismember("Interval HOR", Plot)  == 1
    
    
    %Histogram of Stowood results
    figure 
    hold on
    histogram(Grounds, 10)
    xlabel('Reference subject RR (brpm)')
    ylabel('Estimates') 
    
    % Calculating percentiles 
    Upper = prctile(Grounds,99.5);
    Lower = prctile(Grounds,0.5);
    
    % Plotting percentiles   
    plot([Upper, Upper], [0, 2500], '--k', 'LineWidth', 2)
    plot([Lower, Lower], [0, 2500], ':k', 'LineWidth', 2)
    
    legend('Airflow RR estimates', '99.5th Percentile (33.1 brpm)', '0.5th Percentile (5.9 brpm)')
    hold off 
    
end

if ismember("Results", Plot)  == 1
    
    if Type == 'b'
         
    f = figure;
    f.Position = [100 100 850 200];
    hold on 
    plot(Time_vector, Estimate_value_vector, 'r', 'LineWidth', 1.3)
%     plot(Time_vector, Metronome_vector, 'k' ,'LineWidth', 1.1)
    plot(Time_vector, Ground_truth_vector, 'm' ,'LineWidth', 1.3)
    xlabel('Time (s)')
    ylabel('Estimated RR (brpm)')
    legend('CSI RR estimates ','Airflow RR Estimation')
%     , 'Metronome Target')   
    set(gca,'XMinorTick','on','YMinorTick','on')
    axis([0, length(Metronome_vector), 0, max(Metronome_vector)*1.1])
    

  
    % Highlighting Flagged Stowood
    
    for n = [1:length(Flag)]
        
        % only patching if previously flagged
        if Flag(n) == 1
                   
            % Creating rectangular patch
            h = patch([Time_vector(n)-0.5, Time_vector(n)-0.5, ...
                Time_vector(n)+0.5, Time_vector(n)+0.5],...
                [0, max(Ground_truth_vector) + 5 ,  ...
                max(Ground_truth_vector) + 5, 0 ] , 'r', ...
                'FaceAlpha', 0.3 , 'EdgeColor', 'none');
            
            % Turning Legend off           
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
        end 
        
     
      end 
    hold off
    
    
   
    % plotting for HR parts 
    
    elseif Type == 'a'
    
    figure
    hold on 
    plot(Time_vector, Ground_truth_vector, Time_vector, Estimate_value_vector)
    xlabel('Time (s)')
    ylabel('Estimated RR')
    legend('Stowood', 'CSI', 'Metronome Target')
    title('Reference vs CSI vs Metronome')
    
    % creating patches for unused data 
    if movement == false 
               
                   
            % Creating rectangular patch
            h = patch([15, 15, ...
                225, 225],...
                [0, max(Ground_truth_vector) + 5 ,  ...
                max(Ground_truth_vector) + 5, 0 ] , 'g', ...
                'FaceAlpha', 0.2 , 'EdgeColor', 'none');
            
             h = patch([840, 840, ...
                1020, 1020],...
                [0, max(Ground_truth_vector) + 5 ,  ...
                max(Ground_truth_vector) + 5, 0 ] , 'g', ...
                'FaceAlpha', 0.2 , 'EdgeColor', 'none');
            
            % Turning Legend off           
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
    elseif movement == true 
        
        
             % Creating rectangular patch
            h = patch([360, 360, ...
                540, 540],...
                [0, max(Ground_truth_vector) + 5 ,  ...
                max(Ground_truth_vector) + 5, 0 ] , 'g', ...
                'FaceAlpha', 0.2 , 'EdgeColor', 'none');
            
            % Turning Legend off           
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
     
      end 
    hold off
    
    
    end 
    
end 

if ismember("Matching Results", Plot)  == 1
    
    tiledlayout(2,1)
    
    
        nexttile
    
    hold on 
    plot(Experimental_predictions{length(Experimental_predictions)}{Filter_index}{2}(:,1) -59 ...
        , Experimental_predictions{length(Experimental_predictions)}{Filter_index}{2}(:,2),'m','LineWidth', 1.5);  
    
    plot( Experimental_predictions{length(Experimental_predictions)}{Filter_index}{2}(:,1)-59,...
        Experimental_predictions{length(Experimental_predictions)}{Filter_index}{2}(:,Source_index), 'r' , 'LineWidth', 1.5);  
    
    
    xlabel('Time (s)')
    ylabel('Estimated RR (brpm)')  
    legend('Airflow predictions', 'Unmatched CSI predictions')
    set(gca,'XMinorTick','on','YMinorTick','on')
     axis([150, 600, 5, 25])


    nexttile
    hold on 
    plot(Time_vector, Ground_truth_vector, 'm','LineWidth', 1.5);
    plot(Time_vector, Estimate_value_vector,  'r' , 'LineWidth', 1.5);
    hold off 
    xlabel('Time (s)')
    ylabel('Estimated RR (brpm)')    
     legend('Airflow predictions', 'Matched CSI Predictions')
    set(gca,'XMinorTick','on','YMinorTick','on')
     axis([150, 600, 5, 25])


    
end

if ismember("Results Metronome", Plot)  == 1
    
    if Type == 'b'
        
    figure
    hold on 
%     plot(Time_vector, Ground_truth_vector, Time_vector, Estimate_value_vector, Time_vector, Metronome_vector, 'LineWidth', 1.5)
    plot(Time_vector, Metronome_vector, 'k' ,'LineWidth', 1.5)
    plot(Time_vector, Ground_truth_vector, 'm' ,'LineWidth', 1.5)
    xlabel('Time (s)')
    ylabel('Estimated RR (brpm)')
    legend('Airflow RR Estimation', 'Metronome Target')   
    set(gca,'XMinorTick','on','YMinorTick','on')
    axis([0, length(Metronome_vector), 0, max(Metronome_vector)*1.1])
    
    
  
    % Highlighting Flagged Stowood
    
    for n = [1:length(Flag)]
        
        % only patching if previously flagged
        if Flag(n) == 1
                   
            % Creating rectangular patch
            h = patch([Time_vector(n)-0.5, Time_vector(n)-0.5, ...
                Time_vector(n)+0.5, Time_vector(n)+0.5],...
                [0, max(Ground_truth_vector) + 5 ,  ...
                max(Ground_truth_vector) + 5, 0 ] , 'r', ...
                'FaceAlpha', 0.3 , 'EdgeColor', 'none');
            
            % Turning Legend off           
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
        end 
        
     
      end 
    hold off
    
    end 


end 

 