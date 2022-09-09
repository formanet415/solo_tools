function  [vpars, vpers, dt] = load_pas_vdf_fac(date, deltat, resolution, verbose)
%LOAD_FAC_PAS  Transforms PAS into field aligned coordinates
%   Returns arrays containing the paralell velocity and perpendicular
%   velocity as coordinates and the number of protons.
%   _______________________________________________________
%   Input
%   date = datenum time
%   deltat = time to average the VDFs
%   resolution = resolution of the velocity grid
%   verbose = verbose
%   _______________________________________________________
%   Output
%   vpars = paralel velocity coordinates
%   vpers = perpendicular velocity coordinates
%   dt = VDF data projected onto grid as defined by the coordinates

%#ok<*AGROW>
caa_data_paths
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = 1;
end
if ~exist('resolution', 'var') || isempty(resolution)
    resolution = 100;
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
sinphi = sin(acos(cosphi)); % multiplication factor perpendicular to B
vels = ev_to_vel(energies);

vdf = permute(vdf,[2 4 1 3]);
vpar = [];
vper = [];
nums = [];
for i = 1:length(vels)
    vel = vels(i);
    vpar = [vpar reshape(vel*cosphi,length(ep),[])]; % paralell velocity coordiate
    vper = [vper reshape(vel*sinphi,length(ep),[])]; % perpendicular velocity coordinate
    nums = [nums reshape(vdf(i,:,:,:),length(ep),[])]; % number of protons (original vdf data)
end

% out(1,:,:) = vpar; out(2,:,:) = vper; out(3,:,:) = nums;
vpar = double(vpar); vper = double(vper); nums = double(nums);

if verbose == 1
    disp('interpolating data')
end
c = 299792458;
dt = zeros(length(ep), resolution, resolution);
pamin = min(min(vpar)); pamax = max(max(vpar)); % grid limits for interpolation
pemin = min(min(vper)); pemax = max(max(vper));
% pamax = 0.3*c; pamin = -0.1*c; pemax = 0.3*c; pemin = 0;
vpars = linspace(pamin/c, pamax/c, resolution); % making grid coordinates
vpers = linspace(pemin/c, pemax/c, resolution); 
for i = 1:length(ep)
    F = scatteredInterpolant(vpar(i,:)'/c,vper(i,:)'/c,(nums(i,:)'),'natural','none'); 
    dt(i,:,:) = F(repmat(vpars,resolution,1),repmat(vpers',1,resolution));
end
dt(dt<0) = 0;
dt = squeeze(nanmean(dt));

end

