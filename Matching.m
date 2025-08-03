%%
% File    : Matching.m    
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
%  Matching finds the lag that maximises the cross correlation of vectors a
%  and b. The lag is the time, in seconds, that vector b should be  
%  "advanced by" to optimally match the vectors. Uncommenting code under
%  plot section shows the xcorr as a function of lag. 
%
%
% INPUT
% -----
%
% a: Base vector to which b is matched 
% 
% b: Vector which is matched to a.
% 
% FsRes: Frequency of the outputs.
%
% OUTPUT
% -----
% 
% 
% Lag: he time, in seconds, that vector b should be  
%  "advanced by" to optimally match the vectors


function Lag = Matching(a, b, FsRes)
%%% function takes 2 time series vectors, a and b, and calculates the lag
%%% that b must be delayed in order to maximise the cross correlation
%%% between the two


%applying a median filter to the a & b 
a(:,2) = medfilt1(a(:,2) , 10,  'truncate') ;
b(:,2) = medfilt1(b(:,2) , 10,  'truncate') ;

%detrending
a(:,2) = detrend(a(:,2));
b(:,2) = detrend(b(:,2));

%initializing 
Xcorr_store =[];
m = 1;

Lag_start = -20 ;
Lag_end = 20 ;

%lag represents how much b trails a 
for lag = [Lag_start : 1/FsRes : Lag_end] 
    
    %Initializing 
    Xcorr = 0 ;
    N=0;
    Xcorr1 = [];
    Xcorr2 = [];
    
    %finding the means of the two vectors
%     for n = [1 : length(b)] 
%         
%         %ensuring that it is within the range 
%         if (n - lag) > 0 && (n - lag) <= length(b) && n < length(a)
% 
%             N = N+1;
%             %storing vectors for mean calculation
%             Xcorr1(N) = b(n - lag, 2);
%             Xcorr2(N) = a(n , 2);
%             
%         end
%         
%     end

    N = length(b);
    
    for n = [1 : length(b)] 
        
        %ensuring that it is within the range 
        if (n - lag) > 0 && (n - lag) <= length(b) && n < length(a)
            
            Xcorr = Xcorr + (b(n - lag, 2)-mean(b(:,2))) * (a(n , 2)-mean(a(:,2))); 
            
        end
        
    end 

    Xcorr = ((Xcorr)*N)/((N-1)*var(b(:,2))*var(a(:,2)));

    Xcorr_store(m) = Xcorr;
    
    m = m+1 ;   
end 

Xcorr_store;

%% Plotting
% 
% tt = [Lag_start : 1/FsRes : Lag_end];
% 
% [mxcrr, indxx] = max(Xcorr_store);
% 
% v1 = ([0,  tt(indxx)])
% v2 = ([mxcrr,  mxcrr])
% 
% figure
% plot([Lag_start : 1/FsRes : Lag_end] , Xcorr_store, 'LineWidth', 1.5)
% hold on 
% scatter(tt(indxx), mxcrr, 'r', 'LineWidth', 5);
% plot(v1, v2, 'k--','LineWidth', 2 )
% plot([tt(indxx), tt(indxx)], [mxcrr,  0], 'k--','LineWidth', 2 )
% xlabel('Lag (s)') 
% ylabel('Cross-Correlation')
% ax = gca;
% ax.XAxisLocation = 'origin';
% ax.YAxisLocation = 'origin';
% set(gca,'XMinorTick','on','YMinorTick','on')

%% 
%finding lag that maximises the Xcorr 
[~ , Lag_max ] =  max(Xcorr_store);

%normalizing
Lag_max = Lag_start + Lag_max*FsRes;

%taking negative so lag is should be added to maximise 
Lag =  -1 * Lag_max;
end