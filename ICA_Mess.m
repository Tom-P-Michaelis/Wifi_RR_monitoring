%% Importing data 

%participant number 
Participant = 'P008b';

[CSI, Time_store, ~, ~] = Data_reader_V5(Participant,"PCAP");


 
%% sudo windowing 
CSI = CSI(950*47:980*47,:);

%% splitting into phase and Magnitude 

CSI= abs(CSI); 

%% ICA 

[RR, RR_SNR, ICA_Peak,] = ICA_Filter(CSI)