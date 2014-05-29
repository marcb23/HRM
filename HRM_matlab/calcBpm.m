function [bpm,bpm_str,buffer_pulses,stamps_tail] = calcBpm(stamps_head,...
                        stamps_tail,time_stamps,buffer_pulses,bpm,pos)

% Calculate beats per minute from the last 'buffer' seconds of data. 
% If there is not enough new data keep the last calculated heart rate (or 0,
% the sentinel). If the calculated heart rate is too extreme set it to the
% last legitimate rate (or 0)
    buffer = 10;
    if pos>1 && bpm(pos-1) > 0 
        bpm_str = int2str(bpm(pos-1));
    else
%         disp('bpm_str reset up top');
        bpm_str = '...';
    end
    if(stamps_head > 0 && time_stamps(stamps_head) - ...
                        time_stamps(stamps_tail) > buffer)
        while(time_stamps(stamps_head)-time_stamps(stamps_tail) > buffer)
            buffer_pulses = buffer_pulses - 1;
            stamps_tail = stamps_tail + 1;
        end
        diff = time_stamps(stamps_head) - time_stamps(stamps_tail);
        
        % This method uses less averaging
%         bpm(pos) = 60*((buffer_pulses-1)/diff);

        % alternatively, use the calculated bpm to either increase or
        % decrease the bpm by one
        calc_bpm = 60*((buffer_pulses-1)/diff);
        
        % throw out the calculated bpm if its too ridiculous
        if(calc_bpm > 250 || calc_bpm < 45)
            bpm_str = '...';
            if(pos > 1)
                bpm(pos) = bpm(pos-1);
            else
                bpm(pos) = 0;
            end
        elseif(pos > 1 && bpm(pos-1) > 0)
            if(calc_bpm < bpm(pos-1))
                bpm(pos) = bpm(pos-1)-1;
            elseif(calc_bpm > bpm(pos-1))
                bpm(pos) = bpm(pos-1)+1;
            else
                bpm(pos) = bpm(pos-1);
            end
            bpm_str = int2str(bpm(pos));
        else
            bpm(pos)= calc_bpm;
            bpm_str = int2str(bpm(pos));
        end
        
    elseif pos > 1
        bpm(pos) = bpm(pos-1);
    else
        bpm(pos) = 0;
    end
bpm_str = strcat('HRM Output: BPM = ',bpm_str);
