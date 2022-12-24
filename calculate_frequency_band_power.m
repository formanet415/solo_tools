function calculate_frequency_band_power(fmin,fmax,verbose)
%calculate_frequency_band_power(cdf, f0 [Hz],f1 [Hz]) 
%Calculates the power from TDS data in the range from f0 to f1. 
    
    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 0;
    end
    
    epoch_low_samp = [];
    epoch_high_samp = [];
    band_power_ls = [];
    band_power_hs = [];

    % start browsing tswf data from the beginning of the mission
    t0 = datetime(2020,1,1);
    samemonth = month(t0);
    while t0<datetime(date)
        
        % skipping days with VGAM, EGAM or a type III event
        if gam_detector(t0, verbose) || typeIII_detector(t0, verbose) 
            if verbose==1
                disp('skipping day')
            end
            t0 = t0+1;
            continue
        end
        
        
        % making monthly plots
        if month(t0)~=samemonth
            disp('month done')
            %todo plot the past month of data
        end
        samemonth = month(t0);
        
        
        %loading cdf file
        cdf = tdscdf_load_l2_surv_rswf(t0, 1);
        if isempty(cdf)
            t0 = t0+1;
            continue
        end
        
        
        % calculating the powerspectrum
        for i = 1:length(cdf.epoch)
            uu = convert_to_SRF(cdf,i); % gets rid of nans
            
            % not sure how to do the power spectrum
            % [sp, fq, nav] = make_spectrum(uu(1,:),length(uu),1/cdf.samp_rate(i),f1); %this does something
            try
                p = bandpower(uu',cdf.samp_rate(i),[fmin,fmax]);
            catch exception
              	fprintf('Caught error on index %i, %s \n', i, exception.message)
                continue
            end
            
            if cdf.samp_rate(i)>5e5
                epoch_high_samp(end+1) = cdf.epoch(i); %#ok<AGROW>
                band_power_hs(1:2,end+1) = p; %#ok<AGROW>
            else
                epoch_low_samp(end+1) = cdf.epoch(i); %#ok<AGROW>
                band_power_ls(:,end+1) = p; %#ok<AGROW>
            end
            
        end
        
    
        t0 = t0+1;
        save('temp_backup.mat', 'band_power_hs', 'band_power_ls', 'epoch_high_samp', 'epoch_low_samp')
    end

    save(sprintf('%i_%i_kHz_band_power.mat', fmin, fmax), 'band_power_hs', 'band_power_ls', 'epoch_high_samp', 'epoch_low_samp')

end

