%%
% File    : Linreg.m    
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
% DESCRIPTON:
% ----------
%
% Linreg.m computes the estimated BrPM from the FFT of the windowed signal.
% The peak of the signal in frequency domain is first calculated. The peak
% and adjacent bins are then kept, but all other values are set to 0. This
% leaves a FFT with a single dominant frequency. The inverse of this FFT is
% calulated, leaving a sinosioidal time series function, to which a sin
% wave is fit to via a non linear model. The frequency of the fitted
% function is used to calculate the predicted BrPM.
%
% INPUT
% -----
%
% FFT: The FFT of CSI window. 
% 
% Fs: Sampling frequency. 
% 
% Bin: Number of FFT values adjacent to the peak that should be kept.
%
% OUTPUT
% -----
% 
% omega: Predicted BrPM 

function omega = Linreg(FFT,Fs,Bin)

Y = FFT;
n = length(Y);

FFT_half= (Y(1:floor(n/2)));
Callibrate = (0:1:length(FFT_half)-1)*Fs/(n+1); 

%% Filtering FFT so that only contains RR bins 

FFT_RR = FFT_half;



% lower bound
Cut1 = 5/60; 

% upper bound 
Cut2 = 40/60; 

%finding callibrate index that is greter or equal to cut1 
nn1 = 0;
c = 0;
while c < Cut1
 nn1 = nn1 + 1 ;
 c = Callibrate(nn1);   
end 

FFT_RR([1:nn1-1]) = zeros(nn1-1,1);

%finding callibrate index that is greter or equal to cut1 
nn2 = 0;
c = 0;
while c < Cut2
 nn2 = nn2 + 1 ;
 c = Callibrate(nn2);   
end 

FFT_RR([nn2:length(FFT_RR)]) = zeros( length(FFT_RR) - nn2 +1  ,1);       
        
FFT_half =  FFT_RR ; 

% plotting filtered FFT
% plot(Callibrate, FFT_half)
% title("Filtered FFT")
% 
% pause 

%% peak detection
[~, freq] = max(FFT_half);

%frequency is target frequency 
Frequency = Callibrate(freq);

% %% filter around peak 
% Filter = zeros(length(Y),1);
% 
% %filter 
% for peaks = [freq length(Y)-freq +1]
%    
%   Bin_start = peaks-Bin 
%   Bin_end = peaks+Bin
%   
%   %ensuring binning does not overrun the index
%   if Bin_start <= 0 
%       Bin_start = 1;
%       Bin_end = 1;
%   end 
%   
%   if Bin_end > length(Y) 
%       Bin_end = length(Y);
%       Bin_start = length(Y);
%   end 
%   
% %   %displays 
% %   Bin_start = Bin_start;
% %   Bin_end = Bin_end;
% %   Length_Filter = length(Filter);
%   
%   Filter(Bin_start:Bin_end) = Y(Bin_start:Bin_end);
% 
% end 

%% filter around peak 

peaks = freq;

Filter = zeros(length(FFT_half),1);

%filter   
  Bin_start = peaks-Bin ;
  Bin_end = peaks+Bin;
  
  %ensuring binning does not overrun the index
  if Bin_start < nn1 
      Bin_start = nn1;
      Bin_end =  peaks + (peaks - nn1) ;
  end 
  
  if Bin_start <= 0
      Bin_start = 1;
      Bin_end =  peaks + (peaks - 1) ;
  end 
  
  if Bin_end > nn2
      Bin_end = nn2;
      Bin_start = peaks + (peaks - nn2);
  end
  
%   %displays 
    Bin_start = Bin_start;
    Bin_end = Bin_start;
%   Length_Filter = length(Filter);
  
  Filter(Bin_start:Bin_end) = FFT_half(Bin_start:Bin_end);
  
  
  %recreating FFT 
  
  Y(1:length(Filter)) = Filter ; 
  Y(length(Y) - length(Filter) + 1:length(Y)) = -flip(Filter) ;
  
  %weird divet correction 
  Y(710) = 0 ;

%% plotting 

% delete(gca)
% plot([1:length(Y)], Y)
% title("Binned FFT")
% pause
% 
% % Plotting the filtered FFT 
% delete(gca)
% plot(Callibrate, Filter)
% title("Binned FFT half")
% pause


%% inversing fft 
X_dash = real(ifft(Y));

%% plotting 


%% linreg

%creating model
model = @(b,x)b(1).*sin(b(2).*x+b(3));

%initialising parameters 
Mag = max(X_dash);
W_normal = Frequency *2*pi/Fs;
Offset = 0;

beta0 = [Mag W_normal Offset];

inp = [1:length(X_dash)];
out = X_dash;

fit = fitnlm(inp, out , model, beta0);

% %plotting 
% a1 = fit.Coefficients{:,1}(1);
% a2 = fit.Coefficients{:,1}(2);
% a3 = fit.Coefficients{:,1}(3);
% 
% Y_dash = a1.*sin(a2.*inp+a3);
% 
% delete(gca)
% plot([1:length(X_dash)], X_dash, [1:length(X_dash)], Y_dash)
% title("fit")
% pause


% finding bpm 
omega = (fit.Coefficients{:,1}(2)/(2*pi))*Fs*60;

%% Filtering OMEGA 

Frequency = Frequency*60;

%error metric
Err = abs((Frequency-omega))/Frequency ; 

%thresholding 
if Err > 0.2
    
    %Frequency is metric if linreg fails 
    omega = Frequency;
    
end 


end