%% Intro
% File    : Data_reader_V5.m    
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
% Depending on the value of Type, Data_reader_V5.m ouptuts the Stowood RR ,
% HR or Wi-Fi CSI data in array format. If a converted file has not already 
% been made from PCAP --> CSV, Data_reader saves the CSV.
%
%
% INPUT
% -----
%
% Num is the unique participant number with the experiment section
% identifier attached e.g. P006a
% 
% Type: A string describing the type of data to be read. Options are
% "PCAP", "Airflow" or "Pulse"
% 
% OUTPUT
% -----
%
% CSI, Time_store, Pulse, Airflow: All are arrays containing the relevent
% data components. 
% ________________________________________________________________________
%
% Improvement Aims from previous versions 
% -----
% 
% > Allow Type to dictate the reading object [DONE]
% 
%

function [CSI, Time_store, Pulse, Airflow] = Data_reader_V5(Num, Type)

%% Initialising Data stores 

CSI = [];
Time_store = [];
Pulse = [];
Airflow = [];

%% Data Nomenclature 

%Defining the root of the data wrt C:
Data_root = "C:\Users\tommi\OneDrive\Documents\4th Year Eng\4YP\DATA\Experimental_March_21";

Participant = Num; 

%% Reading PCAPs 

% If defined by Type 
if Type == "PCAP"
    
%getting filenames
PCAP = strcat(Data_root, "\PCAPS\Experimental\", Participant, ".pcap") ;

%Processed PCAP CSV filenames
PCAP_processed_csi = strcat(Data_root, "\Processed PCAPS\CSI\", Participant, ".csv") ;
PCAP_processed_time = strcat(Data_root, "\Processed PCAPS\Time_store\", Participant, ".csv") ;

%% Checking if CSV has already been read 

% if not using reader function to convert PCAP file 

if isfile(PCAP_processed_csi) && isfile(PCAP_processed_time)
           
else
    
%% Converting data

% Reading original PCAP 
[Time_store, CSI, ~] = Reader(PCAP);

%% Saving data 

%saving CSI 
writematrix(CSI, PCAP_processed_csi)

%saving Time 
writematrix(Time_store, PCAP_processed_time)

end

%importing processed PCAP from their saved file location

CSI = readtable(PCAP_processed_csi);
Time_store = readtable(PCAP_processed_time);

%displaying correct reading
display(strcat(num2str(Participant), " data imported"))

if istable(CSI) == 1
 CSI = table2array(CSI);
end

if istable(Time_store) == 1
 Time_store = table2array(Time_store);
end


end


%% Reading RRs

% If defined by Type 
if Type == "Airflow"

%getting filenames
Airflow = strcat(Data_root, "\Stowood Processed\", Participant, "\airflow.csv");

%importing Stowood Airflow data 
Airflow = readtable(Airflow);

if istable(Airflow) == 1
 Airflow = table2array(Airflow);
end


end


%% Reading HRs 

% If defined by Type 
if Type == "Pulse"
    
%getting filenames 
Pulse = strcat(Data_root, "\Stowood Processed\", Participant, "\pulse.csv");

%importing Stowood pulse data 
Pulse = readtable(Pulse);

if istable(Time_store) == 1
 Time_store = table2array(Time_store);
end

end

end