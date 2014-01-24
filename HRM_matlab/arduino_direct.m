clear all;
close all;
a=arduino('COM12');

pos = 1;
scroll_width = 20;
delay = .1;
sample_time = 25;
time = 0;%[1 2 3 4 5];%zeros(1,num_samples);
data = 0;%[200 500 750 500 100];%zeros(1,num_samples);
deriv = 0;
local_min = 0;
local_max = 1000;

% For calculating BPM
bpm = 0;
bpm_str = '';
status = '';
time_stamps = 0;
stamps_head = 0;
stamps_tail = 1;
buffer_pulses = 0;

threshold = 400;
beats = 0;

%  set up the figure
fig = figure(1);
subplot(3,1,1);
hold on;
plot_raw = plot(time,data,'r','LineWidth',2);
title_all = title('HRM Output','FontSize',15);
hold on;
plot_beats = plot(time,beats,'LineWidth',1 );
% xlabel('Time Elapsed','FontSize',10);
ylabel('Arduino','FontSize',10);
axis([0 scroll_width 0 1000]);
grid on;

% plot the calculated heartrate
subplot(3,1,2);
plot_rate = plot(time,bpm,'g','LineWidth',2);
ylabel('Beats Per Minute','FontSIze',10);
axis([0 scroll_width 0 200]);
grid on;


% plot the derivative of the heart rate
subplot(3,1,3);
plot_deriv = plot(time,deriv,'b','LineWidth',2);
% title('Status: Green','FontSize',15);
xlabel('Time Elapsed','FontSize',10);
ylabel('Rate of change','FontSize',10);
axis([0 scroll_width -5 5]);
grid on;

tic
while (1);
    
    data_in = a.analogRead(2);
    current_time = toc;
    
    %amplify the signal from Arduino
    data(pos) = data_in * 2;
    time(pos) = current_time;
    if(data(pos) > threshold)
        beats(pos) = 900;
        
        % only accept a new beat if it's close enough to the established
        % bpm
        if(pos > 1 && bpm(pos-1) > 0)
            validate = (60/bpm(pos-1))/2; 
        end
        if(pos > 1 && beats(pos-1) == 0 && (bpm(pos-1) == 0 ||...
                current_time-time_stamps(stamps_head) > validate))
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = current_time;
            buffer_pulses = buffer_pulses + 1;
        end
    else
        beats(pos) = 0;
    end
    [local_min, local_max] = scaleAxis(pos,data,time);
    threshold = local_max - (local_max-local_min)/2;
  
    if(~ishandle(fig))
        break
    end
    subplot(3,1,1,'align');
    hold on;
    set(plot_raw,'xdata',time,'ydata',data);
    set(plot_beats, 'xdata', time, 'ydata', beats);
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) local_min local_max]);
    else
        axis([0 scroll_width local_min local_max]);
    end

    [bpm,bpm_str,buffer_pulses,stamps_tail] = calcBpm(stamps_head,...
    stamps_tail, current_time, time_stamps, buffer_pulses, bpm, pos);
    set(title_all, 'String', bpm_str);
    
    % plot the beats per minute over time
    subplot(3,1,2,'align');
    set(plot_rate, 'xdata', time, 'ydata', bpm);
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) 0 200]);
    end
    
    % do things to calculate the rate of change of heartrate
    [deriv,status] = calcDeriv(deriv,bpm,pos,current_time,status);
    subplot(3,1,3,'align');
    hold on;
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) -5 5]);
    end
    set(plot_deriv,'XData',time,'YData',deriv);

    
%     set(title_deriv, 'String', status);
    drawnow;
    
    pos = pos + 1;
    pause(delay);
end
