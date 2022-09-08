function plot_pas_vdf_fac(date,deltat, resolution, verbose)
%PLOT_PAS_VDF_FAC Plots PAS FAC VDF
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = 1;
end
if ~exist('resolution', 'var') || isempty(resolution)
    resolution = 100;
end
[x,y,dt] = load_pas_vdf_fac(date, deltat, resolution,verbose);
imagesc(x,y,log10(dt))
xlabel('v_{||}/c','FontSize',16)
ylabel('v_{\perp}/c''FontSize',16)
cmap = [1 1 1; parula(256)];
colormap(cmap)
set(gca, 'YDir', 'normal')
graph = gcf;
graph.Position = [100 100 800 800];
end

