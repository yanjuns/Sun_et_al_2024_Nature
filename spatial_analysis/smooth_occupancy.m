function smoccu = smooth_occupancy(occupancy_curve)
%this function aims to smooth the animal's occupancy curve/map. Grabed from
%Kiah Hardcastle's code compute_2d_tuning_curve_dt (Butler at al., 2019). 
%Yanjun Sun 9/22/22

%new smoothing (very similar to Stensola et al. 2015 and other Moser papers)
%they use 5 x 5 filter for smoothing 3cm bins; Leutgebs use 5 x 5 filter for
%5 cm bins
H = fspecial('gaussian',[5 5],1);
%smooth spikes and occupancy maps separately, then calculate rate map
% smspikes = imfilter(spike_curve,H);
smoccu = imfilter(occupancy_curve, H);
smoccu(occupancy_curve == 0) = 0;
% tuning_curve = smspikes./smoccu;

end