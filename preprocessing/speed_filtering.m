function frames_pass = speed_filtering(speed,thresh_speed_min)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% if  ~exist('thresh_speed_max','var') || isempty(thresh_speed_max)
%     thresh_speed_max = 20; %cm/sec
% end

points_omit = 10;
idx = [];
for ii = points_omit:length(speed)-points_omit
    if speed(ii) >= thresh_speed_min %&& speed(ii) <= thresh_speed_max
        idx = [idx,ii];
    end
end
frames_pass = idx;


end

