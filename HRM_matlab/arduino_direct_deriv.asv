clear all;
close all;
a=arduino('COM10');

pos = 1;
scroll_width = 20;
delay = .1;
sample_time = 25;
time = 0;%zeros(1,num_samples);
data = 0;%zeros(1,num_samples);
deriv = 0;

% For calculating BPM
time_stamps = 0;
stamps_head = 0;
stamps_tail = 1;
buffer_pulses = 0;
buffer = 10;
bpm = 0;
diff = 0;

threshold = 200;
beats = 0;

%  set up the figure
figure(1);
plot_raw = plot(time,data,'r','LineWidth',2);
hold on;
plot_beats = plot(time,beats,'LineWidth',1 );
title_handle = title('HRM Output','FontSize',15);
xlabel('Time Elapsed','FontSize',15);
ylabel('Arduino','FontSize',15);
axis([0 scroll_width 0 1000]);
grid on;
figure(2);
plot_deriv = plot(time,deriv);
grid on;

tic
while (1);
    data(pos) = a.analogRead(1);
    time(pos) = toc;
    deriv(pos) = 0;
    if(data(pos) > threshold)
        beats(pos) = 900;
        validate = (60/bpm)/2;
        if(pos > 1 && beats(pos-1) == 0 && (bpm == 0 ||...
                toc-time_stamps(stamps_head) > validate))
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = toc;
            buffer_pulses = buffer_pulses + 1;
        end
    else
        beats(pos) = 0;
    end
%     point = [tic value];
%     figure(1);
    try
        set(plot_raw,'XData',time,'YData',data);%(time,data);
        set(plot_beats,'Xdata',time,'YData',beats);
        set(plot_deriv,'Xdata',time,'Ydata',deriv);
    catch
        break
    end
    if(time(pos)-scroll_width > 0)
        axis([time(pos)-scroll_width time(pos) 0 1000]);
    end
    drawnow;
    
%     Calculate beats per minute from the last ten seconds of data
    if(stamps_head > 0 && stamps_head - stamps_tail > 5)
        while(toc-time_stamps(stamps_tail) > buffer)
            buffer_pulses = buffer_pulses - 1;
            stamps_tail = stamps_tail + 1;
        end
        diff = time_stamps(stamps_head) - time_stamps(stamps_tail);
        bpm = 60*(buffer_pulses/diff);
        bpm_str = int2str(bpm);
        if(bpm > 250 || bpm < 50)
            bpm_str = '(N/A)';
        end
        new_title = strcat('HRM Output: BPM = ',bpm_str);
        set(title_handle,'String',new_title);
        drawnow;
    end
    pos = pos +1;
    pause(delay);
end
