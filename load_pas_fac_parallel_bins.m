function [outputArg1,outputArg2] = load_pas_fac_parallel_bins(date, deltat, verbose)
%LOAD_PAS_FAC_PARALLEL_BINS Finds the bin closest to the magnetic field
%   Returns the bin data, timetags, energy bins and bin elevation and
%   azimuth.
caa_data_paths
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = 1;
end

ravg = 3; % minimum averaged mag data points
% load ing mag data or LL02 mag data
[b_ep, b_vec] = caadb_get_solo_mag(date, deltat, 'srf');
if isempty(b_ep)
    [b_ep, b_vec] = caadb_get_solo_mag_LL02(date, deltat, 'srf');
    if isempty(b_ep)
        out = [];
        return
    end
end
sr = mean(1./((b_ep(2:end)-b_ep(1:end-1))*86400));
ravg = max(ravg, sr*60);
b_vec = movmean(b_vec, ravg, 2); % Running average over one minute

M = [1 0 0; 0 -1 0; 0 0 -1];
b_vec = M*b_vec; % Transform into PAS coordinates
bn = b_vec./vecnorm(b_vec); % normalize

% loading PAS
[ep, vdf, elev, energies, azim, ~] = caadb_get_solo_swa_pas_vdf(date, deltat);
elevs = repmat(elev, 1, 11); % rows change
azims = repmat(azim',1, 9)'; % columns change

for i = 1:length(ep)
    [~,j] = min(abs(ep(i)-b_ep));
    bni(:,i) = bn(:,j);
end
[x,y,z] = sph2cart(deg2rad(azims), deg2rad(elevs), 1*ones(9,11));
vs(1,:,:) = x; vs(2,:,:) = y; vs(3,:,:) = z;
cosphi = pagemtimes(bni',vs); % multiplication factor paralell to B
phi = acos(cosphi); % angle between B and where the pixels are facing
[phisx, tmp] = max(phi,[],2);
[phisx, azimidx] = max(phisx,[],3);
for i = 1:length(ep)
elevidx(i) = tmp(i,1,azimidx(i));
closestbin(i,:) = vdf(elevidx(i),:,azimidx(i),i);
end




save('bins.mat','energies','closestbin','elevidx','azimidx','ep','phisx')


end