function [beats,time_stamps,stamps_head,buffer_pulses] = beatFinder(...
    beats,pos,bpm,time_stamps,stamps_head,current_time,buffer_pulses)

    % allow the heartbeat to vary by (at most) this amount each time a beat
    % is found
    variation = 1.4;
    diff = 0;
    upper_bound = 250;
    
    % only accept a new beat if it's close enough to the established
    % bpm. EX. when you divide by 1.75 this beat could have come in at
    % <= 1.75 times the rate of the previously established heart rate
    if(pos > 1 && bpm(pos-1) > 0)
        upper_bound = bpm(pos-1)*variation;
        missed_1 = bpm(pos-1)/variation;
        missed_2 = bpm(pos-1)/(2 * variation);
        lower_bound = bpm(pos-1)/(3*variation);
    end
        

    % detect whether this pulse is a new heartbeat
    inst_bpm = 0;
    if(stamps_head > 0)
        diff = current_time-time_stamps(stamps_head);
        inst_bpm = 60/diff;
    end
    if(pos > 1 && beats(pos-1) == 0 && (bpm(pos-1) == 0 ||...
            inst_bpm < upper_bound))
        % detect whether the HRM missed a beat. If it did, add a beat
        % in between the last two
        
        if(bpm(pos-1) > 0 && inst_bpm < missed_1 && ...%validate_upper1 && ...
            inst_bpm > missed_2)
            msg = strcat('add 1, instantaneous: ',int2str(inst_bpm),', last: ',...
                int2str(bpm(pos-1)));
            disp(msg);

            add_time = diff/2 + time_stamps(stamps_head);
            stamps_head = stamps_head+1;
            buffer_pulses = buffer_pulses+1;
            
        % add in 2 beats if it looks like 2 were missed
        elseif(bpm(pos-1) > 0 && inst_bpm < missed_2...% diff > validate_upper2 && ...
            && inst_bpm > lower_bound)
            msg = strcat('add 2, instantaneous: ',int2str(inst_bpm),', last: ',...
                int2str(bpm(pos-1)));
            disp(msg);
            
            add_time1 = diff/3 + time_stamps(stamps_head);
            add_time2 = 2*diff/3 + time_stamps(stamps_head);
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = add_time1;
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = add_time2;
            buffer_pulses = buffer_pulses+2;
        end
        % add a time stamp for the current beat
        stamps_head = stamps_head+1;
        time_stamps(stamps_head) = current_time;
        buffer_pulses = buffer_pulses+1;
    else
        % go here if the most recent beat came in too quickly to be
        % realistic
        if(pos > 1 && beats(pos-1) == 0)
            msg = strcat('ignoring beat, instantaneous: ',int2str(inst_bpm),', last: ',...
                int2str(bpm(pos-1)));
            disp(msg);
        end
    end
end