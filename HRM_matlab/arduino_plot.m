% Phil R
 
% 13/5/2012
 
 
% When a text file is loaded the spaces are removed. Originally I had temperature and humidity logging
 
% which was printed to the serial port like this "H: 63 % T: 22 *C". When its loaded in to freemat
 
% it appears like this "H:63%T:22*C" so we will look for the "%T:" to indicate that the next two digits
 
% are the logged temperature. These will then be extracted and plotted.
 
 
% User Settings
 
clear all %start with fresh workspace
dataPath = 'Z:\Documents\BCH\HRM\hrm_test\';
textFilename = 'arduinolog2.txt';
widthOfTempdata = 3; %eg arduino prints T: 22*C, we are interested in the 22 so that is 2 digits
temperatureDelimiter = 'S'; %This appears before each temperature reading
 
% % Open the file and load it in to
fp = fopen([dataPath textFilename]);
loadedText = fscanf(fp,'%s');
 
% Find Indicies Where TemperatureDelimiter Exists In Loaded String
matchIndex = strfind(loadedText, temperatureDelimiter);
 
% Loop Through and Store Found Temperatue Readings
for i = 1:length(matchIndex)
    loggedData(i) = str2num(loadedText((matchIndex(i)+length(temperatureDelimiter))...
    :(matchIndex(i)+length(temperatureDelimiter))+widthOfTempdata-1));
end
 
% Loop Through and Store Found Humidity Readings
% humidityDelimiter = 'H:'; %This appears before each temperature reading
% matchIndex = strfind(loadedText, humidityDelimiter);
% for i = 1:length(matchIndex)
%     loggedHumidity(i) = str2num(loadedText((matchIndex(i)+length(humidityDelimiter))...
%     :(matchIndex(i)+length(humidityDelimiter))+widthOfTempdata-1));
% end
 
% Plot The Resulting Temperatures
figure(1);clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%random test stuff
rate = 1/length(loggedData);
T = 0.1;
duration = 20;
step = (1/length(loggedData)*20);
t = 0:step:duration-step;
%skip the first 50 samples
start = step*50;

x = cos(2*pi*t*50) + cos(2*pi*t);
plot(t,loggedData);
hold on;
% denominator is the filter time constant
a = T*5;
hpf = filter([1-a a-1], [1 a-1], loggedData);
lpf = filter(a, [1 a-1], loggedData);
xcomp = zeros(1,length(loggedData));
for i=1:length(loggedData)
    if loggedData(i) > 20
        xcomp(i) = 500;
    else
        xcomp(i) = 0;
    end
end
%bpf = 
plot(t,xcomp,'LineWidth',2);
%axis([0 duration*length(loggedData) -2 2]);
% subplot(2,1,1)
%plot(fft(loggedData))
title('HRM Output','fontsize',13,'fontweight','bold');
xlabel('Time (seconds)','fontsize',10,'fontweight','bold');
ylabel('Amplitude','fontsize',10,'fontweight','bold');
axis([start duration 0 999])
 
% Plot The Resulting Humidity
% subplot(2,1,2)
% plot(loggedHumidity,'color','r')
% title('Room Humidity','fontsize',13,'fontweight','bold');
% xlabel('Record Number','fontsize',10,'fontweight','bold');
% ylabel('Temperature [DegC]','fontsize',10,'fontweight','bold');
% axis([0 250 60 70])