%% Importing data 

%participant number 
Participant = 'P011b';

[CSI, Time_store, ~, ~] = Data_reader_V5(Participant,"PCAP");

addpath('C:\Users\tommi\OneDrive\Documents\4th Year Eng\4YP\Matlab\Mauro_Files')


 
%% sudo windowing 
CSI = abs(CSI(730*48:760*48,:));