%participant number 
Participan = 'P011b';

[CSI, Time_store, ~, ~] = Data_reader_V5(Participan,"PCAP");

%% Plotting

x = real(CSI([50:55],20))  ;        
y = imag(CSI([50:55],20)) ;                  
dx = diff(x);
dy = diff(y);
figure
plot(x,y,'ro');
hold on
quiver(x(1:end-1),y(1:end-1),dx,dy, 'b')

%%

figure

m = 20;
my = 20;
x = [real(CSI([58+m:62++my],20)) ; [490 ; 100; -300]]  ;        
y = [imag(CSI([58+m:62+my],20)); [1350 ; 1280; 1160]]  ;   

h=scatter(x,y, 'LineWidth', 4);
axis([min(x)*1.5, max(x)*1.5, min(y)*1.5, max(y)*1.5])

for n = [1:length(y)-1]
    
X= [x(n) x(n+1)];
Y= [y(n) y(n+1)] ;


Pos=get(gca,'position');
xlimits=get(gca,'xlim');
xmin=xlimits(1);
xmax=xlimits(2);
AxesRangeX=diff(xlimits);
FigPosX=Pos(1)+(X-xmin)*(Pos(3)/AxesRangeX);
ylimits=get(gca,'ylim');
ymin=ylimits(1);
ymax=ylimits(2);
AxesRangeY=diff(ylimits);
FigPosY=Pos(2)+(Y-ymin)*(Pos(4)/AxesRangeY);
dX=diff(FigPosX);
dY=diff(FigPosY);
a = annotation('arrow',FigPosX,FigPosY,'headstyle','none', 'Linewidth', 0.1)
a = annotation('arrow',FigPosX-[0 dX./2],FigPosY-[0 dY./2], 'Linewidth', 0.1)
b = a.Color
a.Color = 'red'
end 

ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
% set(gca,'xticklabel',{[]})
% set(gca,'yticklabel',{[]})

xlabel('Real CSI')
ylabel('Imaginary CSI')


%% Randomness determination 

close all 

clear length
len = 30;

data = CSI([50*735:50*735+len*50],28);

Phase = angle(data);
Mag = abs(data);

% for n = [1:length(Mag)] 
%     
%     if Mag(n) >= 1000
%         Mag(n) =  759;
%         
%     end 
% end 

% plotting 
figure 
plot([1:length(Phase)]./50, Phase)
xlabel('Time (s)')
ylabel('CSI Phase (rad)')
axis([0, 30, -1.1*pi, 1.1*pi])


figure
plot([1:length(Mag)]./50, Mag, 'r')
xlabel('Time (s)')
ylabel('CSI Magnitude')
axis([0, 30, min(Mag)*0.99, max(Mag) * 1.01])


%calculating xcorr 
N = 0;
X1 = Mag;
X1 = X1 - mean(X1);

auto = [];

 for lag = [-100:100]
    
     N = 0 ;
     auto1 = 0 
     
     for n = [1:length(X1)]
         
     if n + lag > 0 && n + lag <= length(X1)
         
         auto1 = auto1 + X1(n+lag) * X1(n);
         
         N = N + 1 ;    
     end 
     
       auto(lag + 101) = auto1/N;
       
     end 
         
   
 end 
 
 figure 
 plot(([1:length(auto)]-101)./50, auto, 'r','LineWidth', 2)
 xlabel('Lag (s)')
 ylabel('Autocorrelation')
 axis([-100/50, 100/50, 0, max(auto)*1.1]) 
 
 ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
 
% [autocor,lags] = xcorr(Mag,Mag);
% 
%  plot([1:length(autocor)],autocor)
 
% xlabel('Lag (days)')
% ylabel('Autocorrelation')
% axis([-21 21 -0.4 1.1])


%% 


% Choosing which subjects to analyse. ["all", "b"] selects respiratory data
% from each participant.

 Subjects = ["all", "b"];


% If all subjects are chosen - return string array with every participants
% data

% dummy variable subs 
  Subs = [];    
  

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

for Participant = [1:15]
[CSI, Time_store, ~, ~] = Data_reader_V5(Subjects(Participant),"PCAP");
end 

%% CHOOSE

load('SQI_fusion')

for n = [1:length(SQI_fusion(1,:))]
   Choose(n) =  mean(SQI_fusion(:,n));

end 

figure
plot([1:length(Choose)], Choose, 'LineWidth', 1.5)

deri = strcat('Average CQI of the ', ' n^{' , 'th', '}', ' ranked component') 

deri2 = strcat( '5', '^{' , 'th', '}', ' ranked component with an average CQI of 0.057') 
hold on 

scatter(5, Choose(5), 'r*', 'LineWidth', 2) 

plot([5,5], [0, Choose(5)], 'k--', 'LineWidth', 1.1)
plot([0,5], [Choose(5), Choose(5)], 'k--', 'LineWidth', 1.1)

legend('CQI distribution', deri2)

axis([0, 15, 0, 0.7])

ylabel(deri)
xlabel('Ranked component number')






