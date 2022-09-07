function  dt = load_pas_vdf_fac(date, deltat)
%LOAD_FAC_PAS  Transforms PAS into field aligned coordinates
%   Detailed explanation goes here

ravg = 3; % minimum averaged mag data points
% load ing mag data or LL02 mag data
[b_ep, b_vec] = caadb_get_solo_mag(date, deltat, 'srf');
if isempty(b_ep)
    [b_ep, b_vec] = caadb_get_solo_mag_LL02(date, deltat, 'srf');
    if isempty(b_ep)
        dt = [];
        return
    end
end
sr = mean(1./((b_ep(2:end)-b_ep(1:end-1))*86400));
ravg = max(ravg, sr*60);
b_vec = movmean(b_vec, ravg, 2); % Running average over one minute

M = [1 0 0; 0 -1 0; 0 0 -1];
b_vec = M*b_vec; % Transform into PAS coordinates

% loading PAS
[ep, vdf, elev, energies, azim, extras] = caadb_get_solo_swa_pas_vdf(date, deltat);
elevs = repmat(elev, 1, 11); % rows change
azims = repmat(azim',1, 9)'; % columns change

[x,y,z] = sph2cart(deg2rad(azims),deg2rad(elevs),1*ones(9,11));


end

