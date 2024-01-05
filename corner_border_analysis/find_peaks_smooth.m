%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%find peaks in one dimension with a smoothness factor of 'smoothing_steps'.
%%%If smoothing_steps=1 this function is equivalent to findpeaks.
%%% This version does not consider curve as a circular quantity and may
%%% identify the borders as peaks.
%%% circular is 1 if the curves stand on a circle and ends should be
%%% consider as continuous
function [peaks,values] = find_peaks_smooth(curve,smoothing_steps,circular)

if circular==1
    res=[curve(end-smoothing_steps+1:end), curve, curve(1:smoothing_steps)];
else
    res=[ones(1, smoothing_steps)*(min(curve)-1), curve, ones(1, smoothing_steps)*(min(curve)-1)];
end
raux=ones(1,length(curve));
for i=1:smoothing_steps
    posp=res(smoothing_steps+1+i:end-smoothing_steps+i);
    posm=res(smoothing_steps+1-i:end-smoothing_steps-i);
    raux=raux.*(posp<curve).*(posm<curve);
end
peaks=find(raux);
values=curve(peaks);