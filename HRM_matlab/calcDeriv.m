function [deriv,status] = calcDeriv(deriv,bpm,pos,current_time,status)

%find the derivative of the heartrate as calculated in a parent script

deriv = gradient(bpm);

end