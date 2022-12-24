function out = typeIII_detector(date, verbose)
% Derived from gam_detector
% date - datenum of the beginning of the day (the script will floor it)

% This function returns true if there is a typeIII event from the event
% list happening on this date. 
% Note that the event list only contains
% events with significant EPD electron flux increase and does not contain
% all type III radio emissions.
    
    date = floor(datenum(date));
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 0;
    end
    out = false; 

% Loading data from event list
    load('t3_in_situ_events_V02.mat')
    times = events.rtt(events.rtt>date);

    [~, index] = min(abs(times-date)); % checking if date is close to any Type III events
    if abs(times(index)-date)<1
        out = true;
        if verbose==1
            disp('Date is near to a Type III event')
        end
    end
    
end