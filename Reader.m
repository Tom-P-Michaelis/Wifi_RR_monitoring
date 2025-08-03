function [Time_store, CSI, Time_start] = Reader(filename)
% Takes pcap and converts into a P*T  matrix 

CHIP = '43455c0';          % wifi chip (possible values 4339, 4358, 43455c0, 4366c0)
BW = 80;                % bandwidth
FILE = filename;  % capture file
NPKTS_MAX = 60000000;       % max number of UDPs to process

HOFFSET = 16;           % header offset
NFFT = BW*3.2;          % fft size
p = readpcap();
p.open(FILE);
n = min(length(p.all()),NPKTS_MAX);
p.from_start();
csi_buff = complex(zeros(n,NFFT),0);
Time_store = zeros(n,1);
Time_diff = zeros(n,1);
k = 1;

while (k <= n)
    f = p.next();
    if isempty(f)
        disp('no more frames');
        break;
    end
    if f.header.orig_len-(HOFFSET-1)*4 ~= NFFT*4
        disp('skipped frame with incorrect size');
        continue;
    end
    payload = f.payload;
    H = payload(HOFFSET:HOFFSET+NFFT-1);
    if (strcmp(CHIP,'4339') || strcmp(CHIP,'43455c0'))
        Hout = typecast(H, 'int16');
    else
        disp('invalid CHIP');
        break;
        
    end
    Hout = reshape(Hout,2,[]).';
    cmplx = double(Hout(1:NFFT,1))+1j*double(Hout(1:NFFT,2));
    csi_buff(k,:) = cmplx.';
    
   
    
    time = double(f.header.ts_sec());
    utime = double(f.header.ts_usec());
    Time_stamp = (time + utime/(10^6));
    
    
    if k==1 
       Time_start =  Time_stamp;
    end 
    
    Time_store(k) = Time_stamp - Time_start;
    
    if k~=1 
       Time_diff(k-1) = Time_store(k)- Time_store(k-1);
    end
    
    k = k + 1;
    
end


%exports 
Time_store = [Time_store, Time_diff];
CSI =  csi_buff;  
Time_start = Time_start;
end

