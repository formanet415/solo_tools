function  dt = load_pas_vdf_fac(date, deltat)
%LOAD_FAC_PAS  Transforms PAS into field aligned coordinates
%   Detailed explanation goes here

caa_data_paths
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
bn = b_vec./vecnorm(b_vec); % normalize

% loading PAS
[ep, vdf, elev, energies, azim, extras] = caadb_get_solo_swa_pas_vdf(date, deltat);
elevs = repmat(elev, 1, 11); % rows change
azims = repmat(azim',1, 9)'; % columns change

for i = 1:length(ep)
    [~,j] = min(abs(ep(i)-b_ep));
    bni(:,i) = bn(:,j);
end
[x,y,z] = sph2cart(deg2rad(azims), deg2rad(elevs), 1*ones(9,11));
vs(1,:,:) = x; vs(2,:,:) = y; vs(3,:,:) = z;
cosphi = pagemtimes(bni',vs); % multiplication factor paralell to B
sinphi = sin(acos(cosphi)); % multiplication factor orthogonal to B
vels = ev_to_vel(energies);

vdf = permute(vdf,[2 4 1 3]);
vper = [];
vort = [];
nums = [];
for i = 1:length(vels)
    vel = vels(i);
    vper = [vper reshape(vel*cosphi,length(ep),[])];
    vort = [vort reshape(vel*sinphi,length(ep),[])];
    nums = [nums reshape(vdf(i,:,:,:),length(ep),[])];
end


end

