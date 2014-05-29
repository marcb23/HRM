function [local_min, local_max, threshold, range] = scaleAxis(pos, data, time)

% find the local min and local max of the last buffer_time seconds of the
% graph, which will be used to resize the axes appropriately
buffer_time = 2;
finish = pos;
current_time = time(pos);
while(pos > 1 && time(pos) > current_time - buffer_time)
    pos = pos-1;
end
start = pos;
local_min = data(start);
local_max = data(start);

for i=start:finish
    if(data(i) < local_min)
        local_min = data(i);
    end
    if(data(i) > local_max)
        local_max = data(i);
    end
end

local_min = local_min-5;
if(local_min < 0)
    local_min = 0;
end
local_max = local_max+5;
if(local_max > 1000)
    local_max = 1000;
end
threshold = local_max - (local_max-local_min)/2;
range = local_max-local_min;