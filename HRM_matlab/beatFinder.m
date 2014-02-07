function [beats,time_stamps,stamps_head,buffer_pulses] = beatFinder(...
    beats,pos,bpm,time_stamps,stamps_head,current_time,buffer_pulses)
    % allow the heartbeat to vary by (at most) this amount each time a beat
    % is found
    variation = 1.75;
    diff = 0;
    
    % stop matlab from complaining
    buffer_pulses = buffer_pulses;
        
    % only accept a new beat if it's close enough to the established
    % bpm. EX. when you divide by 1.5 this beat could have come in at
    % <= 1.5 the rate of the previously established heart rate
    if(pos > 1 && bpm(pos-1) > 0)
        validate_lower = (60/bpm(pos-1))/variation;
        validate_upper1 = (60/bpm(pos-1))*variation;
        validate_upper2 = validate_upper1 + 60/bpm(pos-1);
        upper_bound = validate_upper2 + 60/bpm(pos-1);
    end

    % detect whether this pulse is a new heartbeat
    if(stamps_head > 0)
        diff = current_time-time_stamps(stamps_head);
    end
    if(pos > 1 && beats(pos-1) == 0 && ((bpm(pos-1) == 0 ||...
            diff > validate_lower)))
        % detect whether the HRM missed a beat. If it did, add a beat
        % in between the last two
        if(bpm(pos-1) > 0 && diff > validate_upper1 && ...
            diff < validate_upper2)
            add_time = diff/2 + time_stamps(stamps_head);
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = add_time;
            msg = strcat('adding a beat at t=',int2str(10*add_time),...
                ', now=',int2str(10*current_time));
            disp(msg);
            buffer_pulses = buffer_pulses+1;
            
        % optional: add in 2 beats if it looks like 2 were missed
        elseif(bpm(pos-1) > 0 && diff > validate_upper2 && ...
            diff < upper_bound)
            add_time1 = diff/3 + time_stamps(stamps_head);
            add_time2 = 2*diff/3 + time_stamps(stamps_head);
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = add_time1;
            stamps_head = stamps_head+1;
            time_stamps(stamps_head) = add_time2;
            msg = strcat('adding 2 beats at t=',int2str(10*add_time1),...
                ', now=',int2str(10*current_time));
            disp(msg);
            buffer_pulses = buffer_pulses+2;
        end
        stamps_head = stamps_head+1;
        time_stamps(stamps_head) = current_time;
        buffer_pulses = buffer_pulses+1;
    else
        if(pos > 1 && beats(pos-1) == 0)
            msg = strcat('ignoring beat at t=',int2str(10*current_time));
            disp(msg);
        end
    end
end