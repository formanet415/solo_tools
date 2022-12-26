function calculate_frequency_band_power(fmin,fmax,plotit,verbose)
%calculate_frequency_band_power(cdf, f0 [Hz],f1 [Hz]) 
%Calculates the power from TDS data in the range from f0 to f1. 
%
% plotit = 0 - no plotting, 1 - makes the data and then makes plot, 2 - data has been generated previously, just loads them and makes plots
    
    if ~exist('fmin', 'var') || isempty(fmin)
        fmin = 20e3;
    end
    if ~exist('fmax', 'var') || isempty(fmax)
        fmax = 100e3;
    end

    if ~exist('verbose', 'var') || isempty(verbose)
        verbose = 1;
    end
    if ~exist('plotit', 'var') || isempty(plotit)
        plotit = 1;
    end
    
    epoch_low_samp = [];
    epoch_high_samp = [];
    band_power_ls = [];
    band_power_hs_1 = [];
    band_power_hs_2 = [];

    % start browsing tswf data from the beginning of the mission
    t0 = datetime(2020,1,1);
    
    while t0<datetime(date) && plotit~=2
        
        % skipping days with VGAM, EGAM or a type III event
        if gam_detector(t0, verbose) || typeIII_detector(t0, verbose) 
            if verbose==1
                disp('skipping day')
            end
            t0 = t0+1;
            continue
        end
        
        
        % making monthly plots is done at the end

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
            % [sp, fq, nav] = make_spectrum(uu(1,:),length(uu),1/cdf.samp_rate(i),f1); 
            % sum(sp(fq>20e3))*(fq(2)-fq(1))                              <-----   this gives the same result as using bandpower
            
            try % on very few occasions this just does not work. i am using the ostrich method to bugfix this
                p = bandpower(uu',cdf.samp_rate(i),[fmin,fmax]);
            catch exception
              	fprintf('Caught error on index %i, %s \n', i, exception.message)
                continue
            end
            
            if cdf.samp_rate(i)>5e5
                epoch_high_samp(end+1) = cdf.epoch(i); %#ok<AGROW>
                band_power_hs_1(1:2,end+1) = p; %#ok<AGROW>
                
                % higher frequency range
                try     % the try catch might not be needed here
                    p = bandpower(uu',cdf.samp_rate(i),[125e3,200e3]);
                catch exception
                    fprintf('Caught error on index %i, %s \n', i, exception.message)
                    continue
                end
                band_power_hs_2(1:2,end+1) = p; %#ok<AGROW>
            else
                epoch_low_samp(end+1) = cdf.epoch(i); %#ok<AGROW>
                band_power_ls(:,end+1) = p; %#ok<AGROW>
            end
            
        end
        
    
        t0 = t0+1;
        save('temp_backup.mat', 'band_power_hs_1', 'band_power_hs_2', 'band_power_ls', 'epoch_high_samp', 'epoch_low_samp')
    end
    if plotit == 0 || plotit == 1
        save(sprintf('%i_%i_kHz_band_power.mat', fix(fmin/1e3), fix(fmax/1e3)), 'band_power_hs_1', 'band_power_hs_2', 'band_power_ls', 'epoch_high_samp', 'epoch_low_samp')
    elseif plotit == 2
        load('20_100_kHz_band_power.mat') %#ok<LOAD>
    end
    
    if plotit == 0
        return
    end

    for yrs = 2020:year(date) % Plotting all months in plot subdirectory
        for months = 1:12
            lo_dat = [];
            hi_dat_1 = [];
            hi_dat_2 = [];
            emp=1;
            [~,i0] = min(abs(epoch_high_samp-datenum(yrs,months,1)));
            [~,i1] = min(abs(epoch_high_samp-datenum(yrs,months,eomday(yrs,months))));
            if i0==i1
                hi_ep = [];
            else
                emp=0;
                
                hi_ep = epoch_high_samp(i0:i1);
                hi_dat_2 = 1e6*band_power_hs_2(i0:i1);
                plot(hi_ep,hi_dat_2,'g-','Displayname','125-200 kHz power - 524 kHz mode');
                hold on
                
                hi_dat_1 = 1e6*band_power_hs_1(i0:i1);
                plot(hi_ep,hi_dat_1,'r-','Displayname','20-100 kHz power - 524 kHz mode');
            end
            
            [~,i0] = min(abs(epoch_low_samp-datenum(yrs,months,1)));
            [~,i1] = min(abs(epoch_low_samp-datenum(yrs,months,eomday(yrs,months))));
            if i0==i1
                lo_ep = [];
            else
                lo_ep = epoch_low_samp(i0:i1);
                lo_dat = 1e6*band_power_ls(i0:i1);
                plot(lo_ep,lo_dat,'b-','Displayname','20-100 kHz power - 262 kHz mode');
                emp=0;
            end
            if emp == 0
                legend()
                xlim([min([lo_ep,hi_ep]),max([lo_ep,hi_ep])])
                ylim([0,min(0.02,max([lo_dat,hi_dat_1,hi_dat_2]))])
                datetick('x',20,'Keeplimits')
                title(sprintf('Band power plot %i-%i.png',yrs,months))
                hold off
                graph = gcf;
                graph.Position = [100 100 1100 700];
                saveas(graph, ['plots' filesep sprintf('band_power_plot_%i-%i.png',yrs,months)])
                close(graph)
            end
            
            
            
            
            
        end
    end
end

