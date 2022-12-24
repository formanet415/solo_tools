function out = gam_detector(date, verbose)
% This function returns true if there is a gravity assist maneuvere
% happening within 24 hours (before or after).
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 0;
    end
    out = false; 
    within = 48/24; % change this to make the detection tolerance bigger/smaller

% Data retrieved from https://solarorbiter.esac.esa.int/where/ on 2022/12/24
    gams = [datenum(2020,12,27,12,37,0), ...            % VGAM1
        datenum(2021,8,9,4,41,0), ...                   % VGAM2
        datenum(2021,11,27,4,30,0), ...                 % EGAM1
        datenum(2022,9,4,1,25,0), ...                   % VGAM3
        datenum(2025,2,18,20,45,0)];                    % VGAM4 (planned)
        % TODO add future GAMs

    [~, index] = min(abs(gams-date)); % checking if date is close to any GAM
    if abs(gams(index)-date)<within/2
        out = true;
        if verbose==1
            disp('Date is near to a GAM')
        end
    end
    
end