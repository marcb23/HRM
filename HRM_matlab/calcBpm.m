function [bpm,bpm_str,buffer_pulses,stamps_tail] = calcBpm(stamps_head,...
                        stamps_tail,now,time_stamps,buffer_pulses,bpm,pos)

%     Calculate beats per minute from the last ten seconds of data
    buffer = 10;
    bpm_str = '...';
    if(stamps_head > 0 && stamps_head - stamps_tail > 5)
        while(now-time_stamps(stamps_tail) > buffer)
            buffer_pulses = buffer_pulses - 1;
            stamps_tail = stamps_tail + 1;
        end
        diff = time_stamps(stamps_head) - time_stamps(stamps_tail);
        bpm(pos) = 60*((buffer_pulses-1)/diff);
        bpm_str = int2str(bpm(pos));
        if(bpm(pos) > 250 || bpm(pos) < 50)
            bpm(pos) = 0;
            bpm_str = '...';
        end
    else
        bpm(pos) = 0;
    end
bpm_str = strcat('HRM Output: BPM = ',bpm_str);
