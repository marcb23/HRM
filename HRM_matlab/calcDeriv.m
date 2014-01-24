function [deriv] = calcDeriv(bpm,time)

%find the derivative of the heartrate as calculated in a parent script

deriv = gradient(bpm,time);

end