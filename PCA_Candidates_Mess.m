%% Importing data 
close all 

%participant number 
Participant = 'P009b';

[CSI, Time_store, ~, ~] = Data_reader_V5(Participant,"PCAP");

addpath('C:\Users\tommi\OneDrive\Documents\4th Year Eng\4YP\Matlab\Mauro_Files')


CSI = abs(CSI(470*48:500*48,:));
CSI = CSI - mean(CSI);
CSI = detrend(CSI);

%%

% [Candidates, Index_Combined, Weights, Index_store] = PCA_Candidates(CSI, 47, 100);

%% Necking 
% 
% [coeff,scores,latent] = pca(CSI.');
% 
% latent_Sum = [];
% sum = 0 ;
% grad = [];
% 
% for n = [1:length(latent)]
%     
%     sum = sum + latent(n);
%        
%     
% latent_Sum(n) = sum;
% 
% end 
% 
% 
% 
% latent_Sum = latent_Sum./sum;
% 
% for n = [2:length(latent)]
%     
%     grad(n) = (latent_Sum(n)-latent_Sum(n-1))*100;
% 
% 
% end 
% 
% gra = inf;
% r = 1;
% grad(1) = grad(2) ; 
% 
% while gra > 0.1
%     
%     gra = grad(r);
%     
%     r = r+1;
%     
% end 
% 
% 
% figure
% h = plot([1:length(latent_Sum)+1], [0 ,latent_Sum], 'Linewidth', 2);
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% hold on 
% 
% target = 0.975;
% m = 1;
% rel = 0;
% 
% while target > rel 
%     rel = latent_Sum(m);
%     m = m+1 ;
%     
% end 
% 
% h = plot([0,r], [latent_Sum(r), latent_Sum(r)], 'k--', 'LineWidth', 1 );
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% h = plot([r,r], [latent_Sum(r), 0], 'k--', 'LineWidth', 1 );
% h.Annotation.LegendInformation.IconDisplayStyle = 'off';
% scatter(r,latent_Sum(r), 'r*', 'LineWidth', 2)
% 
% 
% gra
% latent_Sum(r)-r*gra
% plot([0, 128], [ latent_Sum(r)-r*gra/100, latent_Sum(r) + (128-r)*gra/100], 'r:', 'LineWidth', 2)
% 
% axis([0,64,0.75,1])
% 
% ylabel('Explained variance')
% xlabel('Dimensions')
% 
% set(gca,'XMinorTick','on','YMinorTick','on')
% legend('97.2% of the variance explained with the first 27 candidates', "gradient of explained variance curve")

%% generating PCA graphs

[coeff,scores,latent] = pca(CSI.');

Fs = 48

    
    HW = hann(length(coeff));    
    scores_hann = coeff .* HW ; 

    %fourier transform
    Length = length(coeff);  
    
    Y_hann = fft(scores_hann);

    %remove sum term 
    FFT_scores = abs(Y_hann(1:floor(Length/2),:));
    Callibrate = (0:1:length(FFT_scores)-1)*Fs/(Length+1); 
    

    tiledlayout(1,3)
    nexttile([1,2])
    
for n = [1:5]
        
            % Selecting Data FFT
            
        FFT = FFT_scores(:,n);
        
        score_proj = scores(:,n);
        coeff_proj = normalize(coeff(:,n), 'range');
        
        col = [0.2*n 0.1+0.1*n 1/(1.5*n)];
        
        if n == 4
            col = [0, 0, 0]
        end 
        
        if n == 3
            col = [1, 0, 0]
        end 
        
        hold on 
        plot([1:length(coeff_proj)]./Fs, coeff_proj + 5-n, 'color', col  )
        
        ylabel('PCA component number')
        xlabel('time (s)')
    
        hold off 
        
end 

axis([0,30, 0, 5.2])

yticks([0.5 1.5 2.5 3.5 4.5])
yticklabels({'5','4','3','2','1'})

axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 4     % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(a)',  'fontweight', 'bold')

nexttile([1,1])

maxxx = 0 ;

for n = [1:5]
        
            % Selecting Data FFT
            
        FFT = FFT_scores(:,n);
        
        score_proj = scores(:,n);
        coeff_proj = normalize(coeff(:,n), 'range');
        
        col = [0.2*n 0.1+0.1*n 1/(1.5*n)];
        
                if n == 4
            col = [0, 0, 0];
        end 
        
        if n == 3
            col = [1, 0, 0];
        end 
        
        hold on 
        h = plot(Callibrate, FFT, 'color', col  ); 
        
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        ylabel('Power (a.u)')
        xlabel('Frequency (Hz)')
        
        [maxval, ind] = max(FFT);
        
        if maxval >= maxxx
        maxxx = maxval;
        Maxind = ind;
        end 
    
end 


axis([0, 1.5, 0, maxxx*1.1])

scatter(Callibrate(Maxind), maxxx, 'k*', 'LineWidth', 3)

legend('FFT peak at 0.23Hz (13.8 brpm)')

axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 0.65   % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(b)',  'fontweight', 'bold')


%% bad value 

% 
%     tiledlayout(1,3)
%     nexttile([1,2])
%     n = 4;
%     
% 
%         
%             % Selecting Data FFT
%             
% %         FFT = FFT_scores(:,n);
%         
%         score_proj = scores(:,n);
%         coeff_proj = normalize(coeff(:,n), 'range');
%         
%         
%         
%         hold on 
%         plot([1:length(coeff_proj)]./Fs, coeff_proj ) 
%         
%         ylabel('CSI (a.u)')
%         xlabel('time (s)')
%     
%         hold off 
%         
% axis([0,30,0,1])
% 
% nexttile([1,1])
% 
%         
%             % Selecting Data FFT
%             
% %         FFT = FFT_scores(:,n);
%         
%         score_proj = scores(:,n);
%         coeff_proj = normalize(coeff(:,n), 'range');
%         
%         
%         
%         hold on 
%         plot(Callibrate, FFT, 'LineWidth', 1) 
%         
%         ylabel('Power (a.u)')
%         xlabel('Frequency (Hz)')
%         
%         [maxy, index] = max(FFT);
%         
%         scatter(Callibrate(index), maxy, '*', 'LineWidth', 2)
%         
%         h = patch([0.0833, 0.083 , ...
%                 0.383, 0.383],...
%                 [0, max(FFT_scores(:,n))*1.2 ,  ...
%                 max(FFT_scores(:,n))*1.2, 0 ] , 'g', ...
%                 'FaceAlpha', 0.3 , 'EdgeColor', 'none');
%         
%         axis([0, 1, 0, max(FFT_scores(:,n))*1.2])
%     
%         hold off 
%         


%%

[coeff,scores,latent] = pca(CSI.');

Filter_config.a = 1
Fs = 48

Time = 0

    
    HW = hann(length(coeff));    
    scores_hann = coeff .* HW ; 

    %fourier transform
    Length = length(coeff);  
    
    Y_hann = fft(scores_hann);

    %remove sum term 
    FFT_scores = abs(Y_hann(1:floor(Length/2),:));
    Callibrate = (0:1:length(FFT_scores)-1)*Fs/(Length+1); 
    

    tiledlayout(1,3)
    nexttile([1,2])
    
    n = 1;
      
            % Selecting Data FFT
            
        FFT = FFT_scores(:,n);
        
        score_proj = scores(:,n);
        coeff_proj = normalize(coeff(:,n), 'range');
        
        RR_CSI = coeff_proj;
        
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
    
    % Plotting filter effect 

        
        t = tiledlayout(2,3)
        
        %plotting unfiltered
        nexttile([1,2])
        plot([1:length(RR_CSI)]./Fs, RR_CSI,  'LineWidth', 1.1)
        
                axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 0.15   % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(a)',  'fontweight', 'bold')
        
        xlabel('Time (s)') 
        ylabel('CSI') 
        set(gca,'XMinorTick','on')
        
        axis([0, 30, min(RR_CSI)*1.1,  max(RR_CSI)*1.1])
        
        % plotting fft 
        nexttile
        
        Y = fft(RR_CSI);
        Y = abs(Y);

        Y_Half = (Y(1:floor(length(RR_CSI)/2),:));
        Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(RR_CSI)+1); 

        plot(Callibrate, Y_Half,  'LineWidth', 1.1 )
        
                axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 0.15   % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(b)',  'fontweight', 'bold')
        
        axis([0, 1 , 0 , max(Y_Half)*1.1])
        
        axis([0, 1 , 0 , max(Y_Half)*1.1])
        
        xlabel('Frequency (Hz)') 
        ylabel('Power') 
        set(gca,'XMinorTick','on')
        
        % plotting filtered
        nexttile([1,2])
        plot([1:length(RR_CSI)]./Fs, RR, 'r', 'LineWidth', 1.1)
        
        axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 0.15   % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(c)',  'fontweight', 'bold')
        
        axis([0, 30, min(RR)*1.1,  max(RR)*1.1])
        
        xlabel('Time (s)') 
        ylabel('CSI') 
        set(gca,'XMinorTick','on')
        
        % next FFT 
        
        nexttile
        
        Y = fft(RR);
        Y = abs(Y);

        Y_Half = (Y(1:floor(length(RR)/2),:));
        Callibrate = (0:1:length(Y_Half)-1)*Fs/(length(RR)+1); 

        plot(Callibrate, Y_Half, 'r',  'LineWidth', 1.1 )
        
        axis([0, 1 , 0 , max(Y_Half)*1.1])
        
        xlabel('Frequency (Hz)') 
        ylabel('Power') 
        set(gca,'XMinorTick','on')
        
        axlim = get(gca, 'XLim');                                           % Get ‘XLim’ Vector
aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
                                                % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 0.15   % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(d)',  'fontweight', 'bold')

        pause   
        delete(gca)
        
        %% SP 
        
               %creating time series 
        t = (0:length(RR)-1)' ./ Fs;
        
        % calculating spectral purity, a low SPI is a pure frequency domain
        [sp_index, ~, ~, ~] = pe.sigproc.spectral_purity(    ...
         RR        , ...
        Fs                  , ...
        t                   , ...
        false                 ...
        )
    
    %% Subcarrier plots 
      Participant = 'P001b';     
   [CSI, Time_store, ~, ~] = Data_reader_V5(Participant,"PCAP");
   
   Fs = 48;


CSI = abs(CSI(657*48:672*48,:));     

tiledlayout(2,3)
nexttile([1,3])
for n = [ 1 ,100 , 200]
plot([1:length(CSI(:,n))]./48, CSI(:,n)-mean(CSI(:,n)), 'LineWidth', 0.8)
hold on 
end 


axis([0, 15 -70 70])
legend('Subcarrier 1', 'Subcarrier 151', 'Subcarrier 200')
xlabel('Time (s)')
ylabel('CSI magnitude (a.u.)')


Filter_config.a = 1;
Time = 1;

nexttile([1,3])
for n = [ 130 ,145 , 200]
    
    RR_CSI = CSI(:,n) ;

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
    
    % Plotting filter effect         
plot([1:length(RR)]./48, RR-mean(RR), 'LineWidth', 1)
hold on 
end 
axis([0, 15 -50 50])
legend('Filtered subcarrier 1', 'Filtered subcarrier 151', 'Filtered subcarrier 200')
xlabel('Time (s)')
ylabel('CSI magnitude (a.u.)')

%% HR, SPO2 

  Subjects = ["all", "b"];

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

a = [];

for n = [1:length(Subjects)]
    
    Data_root = "C:\Users\tommi\OneDrive\Documents\4th Year Eng\4YP\DATA\Experimental_March_21";

Participant = Subjects(n);  


    
%getting filenames
PCAP = strcat(Data_root, "\Stowood Processed\", Participant, "\spo2.csv") ;

    
    a = [a; readtable(PCAP)]; 
    
   
    
end 
    
    
a = a(:,2);

a = table2array(a);

for n = [1:length(a)]
    if n < length(a)
    if a(n) == 0 
        a(n) = 60;
    end 
    end
end 

% Numerical results

Av = mean(a) 
ST = std(a)
    
    %Histogram of results
%     figure 
%     H = histogram(a, 30);
%     xlabel('Mean of Reference RR and CSI RR')
%     ylabel('HR Count')
%     
%     J = boxplot(a, 'orientation', 'horizontal', 'Positions', (max(H.Values)*1.2),'Widths',max(H.Values)*0.1, 'Whisker',1 );
    
    
    %Histogram of Stowood results
    
    figure 
    
    H = histogram(a, 30)
    xlabel('HR (bpm)')
    ylabel('HR count')
      
    set(gca,'XMinorTick','on')
    
    hold on     
    
    J = boxplot(a, 'orientation', 'horizontal', 'Positions', (max(H.Values)*1.2),'Widths',max(H.Values)*0.1, 'Whisker',1 );
    
    set(J,{'linew'},{1})
    yticks([])
    
    % resetting y ticks 
    tik = [0:20000:max(H.Values)*1.2];
    yticks(tik)
    
    
    labs = {};
    for n = [1:length(tik)]
        
  
        labs{n} = char(num2str(tik(n)));
        
    end 
    
    yticklabels(labs)
    
    axis([40,120,0, max(H.Values)*1.5]) 
    
    hold off 
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    figure 

    histogram(a, 100)
    xlabel('SpO2 (%)')
    ylabel('SpO2 count')
      
    set(gca,'XMinorTick','on')
    
    hold on     
    
    J = boxplot(a, 'orientation', 'horizontal', 'Positions', (max(H.Values)*1.2),'Widths',max(H.Values)*0.1, 'Whisker',1 );
    
    set(J,{'linew'},{1})
    yticks([])
    
    % resetting y ticks 
    tik = [0:50000:max(H.Values)*1.2];
    yticks(tik)
    
    
    labs = {};
    for n = [1:length(tik)]
        
  
        labs{n} = char(num2str(tik(n)));
        
    end 
    
    yticklabels(labs)
    
    axis([90,100,0, max(H.Values)*0.9]) 
    
    hold off 
    
    


    %% Metronome ladder
    
    
       %recreating the target RR with arbitrary start time 
    Metronome = [6 : 3 : 33];
    
    FsRes = 1;

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
    
    
    plot([1:length(Metronome_target(:,2))]./60, Metronome_target(:,2), 'k', 'LineWidth', 2)
    axis([0,20,0,35])
    xlabel('Time (minutes)')
    ylabel('Metronome frequency (cycles/minute)')


%% Mouth breath

Fs = 48;
Filter_config.a = 1;
Time = 1;

%participant number 
Participant = 'P004b';

clear Air

[~, ~, ~, Air] = Data_reader_V5(Participant,"Airflow");


Air = Air(1062*32:1087*32,2);
Air = Air - mean(Air);
% Air = detrend(Air);

tiledlayout(2,3)
nexttile([1,3])
plot([1:length(Air)]./32, Air, 'LineWidth', 1)
xlabel('Time (s)')
ylabel('Airflow (a.u.)')
axis([0,25, min(Air)*1.1, max(Air)*1.1])

aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
axlim = get(gca, 'XLim');                                                  % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 3.8   % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(a)',  'fontweight', 'bold')


[CSI, Time_store, ~, ~] = Data_reader_V5(Participant,"PCAP");

CSI = abs(CSI(1062*48:1087*48,160));
CSI = CSI - mean(CSI);
CSI = detrend(CSI);

RR_CSI = CSI ;

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
    
    CSI = RR;

nexttile([1,3])
plot([1:length(CSI)]./48, CSI, 'r', 'LineWidth', 1)
xlabel('Time (s)')
ylabel('CSI (a.u.)')
axis([0,25, min(CSI)*1.1, max(CSI)*1.1])

aylim = get(gca, 'YLim');                                           % Get ‘YLim’ Vector
axlim = get(gca, 'XLim');                                                  % Get ‘legend’ 'Position' Vector
x_txt = min(axlim) - 3.8  % Set ‘x’-Coordinate Of Text Object
y_txt = min(aylim) + diff(aylim)*0.5       % Set ‘y’-Coordinate Of Text Object
text(x_txt, y_txt, '(b)',  'fontweight', 'bold')

   %% Subcarrier plots 2
      Participant = 'P001b';     
   [CSI, Time_store, ~, ~] = Data_reader_V5(Participant,"PCAP");
   
   Fs = 48;


CSI = abs(CSI(657*48:672*48,:));     

tiledlayout(3,3)

for n = [ 1 ,100 , 200]
    nexttile([1,3])
    if n == 1 
 plot([1:length(CSI(:,n))]./48, zeros(length(CSI(:,n)),1),  'r', 'LineWidth', 1)       
    elseif n == 100
plot([1:length(CSI(:,n))]./48, CSI(:,n)-mean(CSI(:,n)), 'k', 'LineWidth', 0.8 )
    elseif n == 200
plot([1:length(CSI(:,n))]./48, CSI(:,n)-mean(CSI(:,n)),  'b', 'LineWidth', 0.8)
    end 
axis([0, 15 -70 70])
xlabel('Time (s)')
ylabel('CSI magnitude (a.u.)')
hold on 
end 






