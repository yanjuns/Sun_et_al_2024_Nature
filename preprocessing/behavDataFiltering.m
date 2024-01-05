function frames_pass = behavDataFiltering(behav,speed,thresh_speed_min,thresh_speed_trans,outliers_removal,factor_upper,factor_lower,thresh_edge,points_omit)
%Function determines indices that fulfill the following criteria
%The indices are mode of contiguous samples
%At each index, the speed is greater than thresh_speed_min
%in each contiguous sample, the maximum speed exceeds thresh_speed_trans

% thresh_edge: defining the edge of the track, which is used for eliminating outliers in the middle of the track. default = 0.15 (i.e. 15% tracklength);
% outliers_removal: a boolean that determines whether to keep or remove outliers.default=TRUE (removes the outliers).
% points_omit: the number of data points (frames) that will be omitted in the start and end during the whole movement; default = 20;

if  ~exist('points_omit','var') || isempty(points_omit)
    points_omit = 10;
end
if  ~exist('thresh_edge','var') || isempty(thresh_edge)
    thresh_edge = 0.15;
end
if  ~exist('factor_lower','var') || isempty(factor_lower)
    factor_lower = 1.05;
end
if  ~exist('factor_upper','var') || isempty(factor_upper)
    factor_upper = 1.5;
end
if  ~exist('outliers_removal','var') || isempty(outliers_removal)
    outliers_removal = true;
end

tracklength = behav.trackLength;
xPosition = behav.position(:,1); yPosition = behav.position(:,2);

i=points_omit;
indices=[];
while i <= length(speed)-points_omit
    if abs(speed(i))>thresh_speed_min
        j=i;
        while j<=length(speed)-points_omit & abs(speed(j))>thresh_speed_min
            j=j+1;
        end
        if max(speed(i:j))>thresh_speed_trans
            indices=[indices,i:j];
        end
        i=j;
    else
        i=i+1;
    end
end

largeyindices=[];
midindices=[];
for i=indices
    if xPosition(i)> thresh_edge*tracklength & xPosition(i)<tracklength-thresh_edge*tracklength
        midindices=[midindices,i];
    end
end
% figure
% boxplot(yPosition(midindices))
meanMid=mean(yPosition(midindices));
stdMid=std(yPosition(midindices));
for i=midindices
    if yPosition(i)>meanMid+factor_upper*stdMid || yPosition(i)<meanMid-factor_lower*stdMid
        largeyindices=[largeyindices,i];
    end
end
if outliers_removal == true
    indices=setdiff(indices,largeyindices);
end
frames_pass = indices;