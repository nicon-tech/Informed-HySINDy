function [value, isterminal, direction] = events_FIEH(t,y,stroke,Yi)
% y(2) is velocity
% y(1) is position
value = abs(y(1))-(stroke/2);
isterminal = 1;
direction = 0;

% value(i) is the value of the ith event function.

% isterminal(i) = 1 if the integration is to terminate at a zero of this
% event function. Otherwise, it is 0.

% direction(i) = 0 if all zeros are to be located (the default). A value
% of +1 locates only zeros where the event function is increasing, and -1
% locates only zeros where the event function is decreasing.