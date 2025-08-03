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



function [RR, HR ,RR_SNR , HR_SNR, PCA_Peak] = PCA_Filter(CSI, Fs, Time)
    %%%Function takes all subcarriers from a window of the CSI. Using PCA it
    %%%aims to detect if there are frequency components corresponding to the HR
    %%% or RR and outputs these principle components. 
     
    %% PCA decomposition 
    
    [~, scores] = pca(CSI);

    %% Fourier transform of scores
    
    % Hann window 
    HW = hann(length(scores));
    
    scores_hann = scores .* HW ; 
    

    %fourier transform
    n = length(scores);                         
    Y = fft(scores);
    Y_hann = fft(scores_hann);

    %remove sum term 
    FFT_scores = abs(Y(1:floor(n/2),:));
    FFT_hann = abs(Y_hann(1:floor(n/2),:));
    Callibrate = (0:1:length(FFT_scores)-1)*Fs/(n+1); 
    
    %% Selecting the best PCs from the first 10 components 

    %initialising 
    HR_counter = 1;
    RR_counter = 1;

    RR_PCA = [] ;
    HR_PCA = [];

    SNR_MAX = 0 ;
    Score_max = 0;
    SNR_RR_Max = 0;
    Peak_Prime = 0;

    %% Looping over first PCs
    for n = [1:150] 
     
        FFT = FFT_hann(:,n) ;

        %calculating baseline mean devaition for the first 2.5 Hz   
        
        Cutoff = round((2.5*(length(FFT)-1)/Fs))*2;
        
        MD = mean(FFT(1:Cutoff)) ;
   
%        MD = mean(FFT) ;
        %% Filtering FFT so that only contains RR bins 
%         
%         FFT_RR = FFT;
%         
%         % lower bound
%         Cut1 = 5/60; 
%         
%         % upper bound 
%         Cut2 = 50/60; 
%         
%         %finding callibrate index that is greter or equal to cut1 
%         nn = 0;
%         c = 0;
%         while c < Cut1
%          nn = nn + 1 ;
%          c = Callibrate(nn);   
%         end 
%         
%         FFT_RR([1:nn-1]) = zeros(nn-1,1);
%         
%         %finding callibrate index that is greter or equal to cut1 
%         nn = 0;
%         c = 0;
%         while c < Cut2
%          nn = nn + 1 ;
%          c = Callibrate(nn);   
%         end 
%         
%         FFT_RR([nn:length(FFT_RR)]) = zeros( length(FFT_RR) - nn +1  ,1);       
%         
        FFT_RR = FFT;

        %% Selecting RR PCs 
        
        [Peak, Freq] = max(FFT_RR);

        Frequency = Callibrate(Freq);

        %% Finding Bimodality
        
        [pks, locs] = findpeaks(FFT_RR);
        
        [Big_Peak, Big_loc] = max(FFT_RR);
        
        Peak_Secondary = 0; 
        Pk_max = 0;
        
        for pk = [1:length(pks)]
            
            if pks(pk) >= Peak_Secondary && abs(locs(pk) - Big_loc) >= 5
                
                Peak_Secondary = pks(pk);       
                Location_Secondary = locs(pk);
                Pk_max = pk;
                
            end 
        
        end 
        
        Lag_Pos = Location_Secondary;
        Rival_Xcorr_max = Peak_Secondary;
        Master_Xcorr = Big_Peak;
        
        % rival peak is relevant if 2/3 of the size 
        if Rival_Xcorr_max >= (1/3)* Master_Xcorr && Frequency <= 0.15          
            Distribution = "Bimodal";      
        else 
            Distribution = "Unimodal";
        end 

       %% PCA Selection 
       
        %if the peak is found at a valid RR frequency
       if Frequency > 0.09 && Frequency < 0.6 && Distribution ~= "Bimodal"

                % if peak of FFT is far greater than white noise then it is significant
                SNR_RR = (Peak)/(MD);
                
                
                % Creating a score that is used to calculate the optimal
                % PCA choice 
                
                % Score penalises low varience as well as rewarding high SNR 
                Score = SNR_RR; 

           if  Score > Score_max

                %storing potential RR PCAs 
                RR_PCA = scores(:,n);
                
                Chosen_PCA = n;
                
                Peak_Prime = Frequency;

                Score_max = Score ;
                
                SNR_RR_Max = SNR_RR;
                
            end              
       end 
       
       %% Selecting HR PCs
       
%        %filtering around HR
%        
%        HR_slice = 1+([0.8, 2] ./ (Fs/(length(scores)+1)));
%        
%        HR_FFT = FFT(round(HR_slice(1)):round(HR_slice(2)) , 1 ); 
%        
%        [HR_Peak, ~] = max(HR_FFT);
%        
%        SNR = HR_Peak/MD;
%        
%        if  SNR > SNR_MAX
% 
%                 %storing potential HR PCAs 
%                 HR_PCA = scores(:,n);
% 
%                 SNR_MAX = SNR ;
%        end 
%      
       %% plotting all PCs 
      
%       
%       if Time >= 456 && Time <= 457
%               
%        tiledlayout(2,1)
%        
%        nexttile 
%        plot(Callibrate, FFT_RR)
%        axis([0, Callibrate(length(Callibrate))/10 , 0 , max(FFT_RR)*1.1])
%        title(strcat("PCA Number ", num2str(n), " FFT"))
%        xlabel("Frequency (Hz)")
%        set(gca,'xMinorTick','on')
%        ylabel("FFT Magnitude")
%        set(gca,'yMinorTick','off')
%        
%        %plotting fft peak 
%        hold on
%        [ymax , xmax] = max(FFT_RR);
%        scatter(Callibrate(xmax), ymax)
%        hold off
%        
%        %plotting secondary peak
%        if Distribution == "Bimodal"
%            
%             hold on
%             scatter(Callibrate(locs(Pk_max)), FFT_RR(locs(Pk_max)), 'k')
%             hold off 
%             
%        end 
%            
%        
%        nexttile 
%        plot([1:length(scores(:,n))]./Fs, scores(:,n))
%        title(strcat("PCA Number ", num2str(n), " Time domain"))
%        xlabel("Time (s)")
%        ylabel("CSI Magnitude")
%        
%        end 
%        
%        pause 
  
    end 
      
%     %% Plotting all PCs at once 
%     
%      if Time >= 800 && Peak_Prime*60 <= 10 
%          
%      tiledlayout(5,2)
%        
%      for component = [1:10]
%          
%        FFT_RR = FFT_hann(:,component);
%          
%        nexttile 
%        plot(Callibrate, FFT_RR)
%        axis([0, Callibrate(length(Callibrate))/12 , 0 , max(FFT_RR)*1.1])
%        title(strcat("PCA Number ", num2str(component), " FFT"))
%        xlabel("Frequency (Hz)")
%        set(gca,'xMinorTick','on')
%        set(gca,'yMinorTick','on')
%        ylabel("FFT Magnitude")
%        set(gca,'yMinorTick','off')
%        
%        %plotting fft peak 
%        hold on
%        [ymax , xmax] = max(FFT_RR);
%        scatter(Callibrate(xmax), ymax, 300, '.')
%        hold off   
%        
%        % if best component change colour
%        if component == Chosen_PCA
%            ax = gca;
%            ax.XColor = 'red';
%            ax.YColor = 'red';
%        end      
%      end
%      
%       
%       
%       pause
%      end 
     
     %% Plotting similar CSI 
%      
%      if Time >= 600 && Time <= 601
%          figure
%          hold on
%      for num = [10, 140, 190, 210,33] 
%          
%          % filtering CSI / detrending 
%          
%          CSI2 = CSI([1:round(length(CSI(:,1))/3)] , num);
% 
%         Filt_LP_RR = designfilt('lowpassfir', 'FilterOrder', 100, ...
%             'StopbandFrequency', 1.1, 'PassbandFrequency', .5, 'SampleRate', 50);
% 
%         %double pass lp filter
%         FIR_filt = filtfilt(Filt_LP_RR, CSI2);
%          
%         FIR_filt = detrend(FIR_filt);
%          
%       
%         plot([1:length(FIR_filt)]./Fs, FIR_filt)  
%          
%         xlabel("Time (s)") 
%         ylabel("CSI Magnitude")
%                 
%          
%      end   
%      hold off 
%      pause 
%      end 
%       
    
    
    %% Outputting PCA analysis 

    %if empty use largest PCA
    if isempty(RR_PCA) == 1 
        RR_PCA = scores(:,1);
        Peak_Prime = 0;
    end 

    if isempty(HR_PCA) == 1  
        HR_PCA(:,1) = scores(:,1);
    end 


    RR = RR_PCA;
    HR = HR_PCA(:,1);
    
    RR_SNR = SNR_RR_Max;
  
    HR_SNR = 1;
    

    %% plotting Outputs 
    
    % tiledlayout(2,1)
    % 
    % nexttile
    % plot([1:length(RR_PCA)]./50, RR_PCA)
    % 
    % nexttile
    % plot([1:length(RR_PCA)]./50, HR_PCA)
    % pause 
    
    %% Outputting 

    PCA_Peak = Peak_Prime;

end 