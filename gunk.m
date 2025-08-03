     Fs = 32 ;

     X = [1:1000];
 RR = rand(1,10000) ;
 
RR = sin(X) ;

%creating time series 
        t = (0:length(RR)-1)' ./ Fs;
        
        % calculating spectral purity, a low SPI is a pure frequency domain
        [sp_index, ~, ~, ~] = pe.sigproc.spectral_purity(    ...
         RR        , ...
        Fs                  , ...
        t                   , ...
        false                 ...
        );
    
        % inversing SPI so high SPI is a good signal 

  sp_index = sp_index
  
  %%
  
       X = [1:1000];
      RR = sin(X/15).' ;
 
  

    HW = hann(length(RR));    
    RR_hann = RR .* HW ; 

    
    %fourier transform
    Length = length(RR);  
    
    Y_hann = abs(fft(RR_hann));

    %remove sum term 
    FFT = abs(Y_hann(1:floor(Length/2),:));
    Callibrate = (0:1:length(FFT)-1)*Fs/(Length+1); 
    
    plot(Callibrate, FFT)
  
      Noise  = mean(FFT);
        
        % Taking valid signal range as between 0.08 and 0.60  
        S1 = 0.08;
        S2 = 0.6;
        
        [~, SI1] = min(abs(Callibrate - S1 - 0.04));
        
        % increasing range to account for spilling 
        
        [~, SI2] = min(abs(Callibrate - S2  +  0.2));      
        
        Signal = mean(FFT(SI1:SI2) );
        
        SNR = Signal/Noise 