% clear all variables and open figures, create an arduino object
clear all;
close all;
a=arduino('COM12');

% basic data plotting variables
range = 1000;   % initial range of the y-axis
default_scale = 1;  % used to scale up or down the analogRead
scale_factor = default_scale;
pos = 1;            % index for the data, beats, time, bpm, and deriv arrays
scroll_width = 20;  % width of the graph window
delay = .1;         % delay between each sample from the arduino
time = 0;           % array for storing raw data and sample times
data = 0;
deriv = 0;          % array to store the rate of change of heart rate

local_min = 0;      % these will be scaled up and down in scaleAxis to 
local_max = 1000;   % match the range of values from arduino

% For calculating BPM, mostly used in calcBpm
bpm = 0;            % array of heart rate calculations
buffer_pulses = 0;  % this will be used in calcBpm to average out a heart rate
bpm_str = '';       % current bpm stored as string
status = '';        % not used yet
time_stamps = 0;    % array of time stamps for all registered heart beats
stamps_head = 0;    % position of the most recent relevant beat
stamps_tail = 1;    % position of the least recent relevant beat
threshold = 400;    % threshold for registering a heart beat, gets adjusted
                    % based on local_max and local_min
beats = 0;          % used to plot the beats, contains all 900s and 0s

% used to save a heart rate value for the csv file every second
second_counter = 1;

%  set up the figure
fig = figure(1);
subplot(3,1,1);
hold on;
plot_raw = plot(time,data,'r','LineWidth',2);
title_all = title('HRM Output','FontSize',15);
hold on;
plot_beats = plot(time,beats,'LineWidth',1 );
ylabel('Arduino','FontSize',10);
axis([0 scroll_width 0 1000]);
grid on;

% plot the calculated heartrate
subplot(3,1,2);
plot_rate = plot(time,bpm,'g','LineWidth',2);
ylabel('Beats Per Minute','FontSIze',10);
axis([0 scroll_width 50 150]);
grid on;


% plot the derivative of the heart rate
subplot(3,1,3);
plot_deriv = plot(time,deriv,'b','LineWidth',2);
% title('Status: Green','FontSize',15);
xlabel('Time Elapsed','FontSize',10);
ylabel('Rate of change','FontSize',10);
axis([0 scroll_width -5 5]);
grid on;

% open up osc communnications with the game
u = udp('127.0.0.1',9000);
fopen(u);

% start a timer at 0 and record the time at which the timer starts
tic
start_time = clock;
while (1);
    % run until the figure gets closed
    if(~ishandle(fig))
        disp('figure closed, breaking out of loop');
        break
    end
    % scale up the signal if its too quiet
    % before making this work add a DC adjustment to take care of noise
%     if(range < 200 && threshold < 600)
%         disp('scaling up, range='); disp(range);
%         scale_factor = scale_factor+.1;
%     elseif(range > 500)
%             scale_factor = default_scale;
%     end

% read in data from arduino and store the current time
    data_in = a.analogRead(2);
    current_time = toc;
    
    % stop sampling if arduino stops sending data
    if(isempty(data_in))
        break
    end
    
    % if there is a defined scale factor, use it to scale input up or down
    data(pos) = data_in * scale_factor;
    time(pos) = current_time;
    
    % if a possible beat is detected run beatFinder
    if(data(pos) > threshold)
        beats(pos) = 900;
        
        [beats,time_stamps,stamps_head,buffer_pulses] = beatFinder(beats,pos,...
        bpm,time_stamps,stamps_head,current_time,buffer_pulses);
    else
        beats(pos) = 0;
    end
    
    % scale the axis to the incoming data. Also sets threshold based on the
    % range of recent arduino readings
    [local_min, local_max, threshold, range] = scaleAxis(pos,data,time);
    
    % update the raw data subplot
    subplot(3,1,1,'align');
    hold on;
    set(plot_raw,'xdata',time,'ydata',data);
    set(plot_beats, 'xdata', time, 'ydata', beats);
    
    % if necessary, scroll the graph window
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) local_min local_max]);
    else
        axis([0 scroll_width local_min local_max]);
    end

    % calculate the heartrate
    [bpm,bpm_str,buffer_pulses,stamps_tail] = calcBpm(stamps_head,...
                    stamps_tail, time_stamps, buffer_pulses, bpm, pos);
    set(title_all, 'String', bpm_str);
    
    % send data to the game over osc
    oscsend(u,'/HXM','iiiis',0,bpm(pos),0,0,'black');
    
    % plot the beats per minute over time
    subplot(3,1,2,'align');
    set(plot_rate, 'xdata', time, 'ydata', bpm);
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) 50 150]);
    end
    
    % do things to calculate the rate of change of heartrate
    % eventually something in this section should determine the status
    % (green, yellow, red) from available data
    % plot the rate of change of the heart rate
    [deriv] = calcDeriv(bpm,time);
    subplot(3,1,3,'align');
    hold on;
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) -5 5]);
    end
    set(plot_deriv,'XData',time,'YData',deriv);
%     set(title_deriv, 'String', status);

    drawnow;
    
    % save the bpm with a timestamp to eventually be exported to a csv
    if(current_time > second_counter)
        output(second_counter,1) = current_time;
        output(second_counter,2) = round(bpm(pos));
        second_counter = second_counter+1;
    end
    
    % increment position in the arrays every time
    pos = pos + 1;
    pause(delay);
end

%close osc communications
fclose(u);

% save output to file
start_time = datestr(start_time,'HHMMSS');
filename = strcat('a_out_',start_time,'.csv');
csvwrite(filename,output);