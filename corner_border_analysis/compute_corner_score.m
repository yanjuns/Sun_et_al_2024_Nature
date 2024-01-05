function cmetrics = compute_corner_score(ratemap,env_coor,mask,penalty,threshlevel)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if ~exist('threshlevel','var')||isempty(threshlevel)
    threshlevel = 0.3;
end
if ~exist('penalty', 'var') || isempty(penalty)
    penalty = true;
end
if ~exist('mask', 'var') || isempty(mask)
    mask = false;
end

if isempty(env_coor)
    cmetrics = struct;
    cmetrics.env_coor = [];
    cmetrics.field_coor = [];
    cmetrics.field2corner_dist = [];
    cmetrics.cornerscore_ea = [];
    cmetrics.cornerscorep = [];
    cmetrics.cornerscore = [];
else
    cmetrics = struct;
    cmetrics.env_coor = env_coor;
    cmetrics.field_coor = cell(1,length(ratemap));
    cmetrics.field2corner_dist = cell(1,length(ratemap));
    cmetrics.cornerscore_ea = cell(1,length(ratemap));
    cmetrics.cornerscorep = cell(1,length(ratemap));
    cmetrics.cornerscore = NaN(length(ratemap),1);
    
    for ii = 1:length(ratemap)
        rtmap = ratemap{ii};
        rtmap(isnan(rtmap)) = 0;
        pkfr = max(max(rtmap));
        %filter out noisy signals
        thresh = threshlevel*pkfr; %0.3 for experiment A, 0.4 for experiment C
        rtmap(rtmap < thresh) = 0;
        %slightly smooth the map for accuracy
        rtmap = imgaussfilt(rtmap,0.5);
%         %under convex condition, use masked region to compute
%         if mask
%             xy = env_coor(1:end-1,:);
%             msk = poly2mask(xy(:,1),xy(:,2),size(rtmap,1),size(rtmap,2));
%             msk = imresize(msk, (max(size(msk))+4)/max(size(msk)));
%             msk = msk(3:end-2,3:end-2);
%             rtmap(~msk) = 0;
%         end
        %find the local maxima of the ratemap and smooth it again
        BW = imregionalmax(rtmap);
        %under convex condition, use masked region to compute
        if mask
            xy = env_coor(1:end-1,:);
            msk = poly2mask(xy(:,1),xy(:,2),size(rtmap,1),size(rtmap,2));
            %the folowing two lines was trying to make the mask slightly
            %larger so that will cover more spaces for field detection. 
%             msk = imresize(msk, (max(size(msk))+4)/max(size(msk)));
%             msk = msk(3:end-2,3:end-2);
            BW(~msk) = 0;
        end
        [row,col] = find(BW == 1);
        field_coor = [col,row];
        %calculate pairwise field2corner_dist between fields and corners, the last row is
        %the field2corner_dist between fields and the centroid.
        D = pdist2(env_coor,field_coor);
        %corner score for a given field is defined as the field2corner_dist to the 
        %centriod minus the min field2corner_dist to a corner and devided by the
        %sum of the two. 
        %The value range of a corner score is -1 to 1, with 1 defined a
        %field that is right at a corner. 
        cornerscore_ea = (D(end,:) - min(D(1:end-1,:)))./(D(end,:) + min(D(1:end-1,:)));
        %this is trying to penalize the cornerscore_ea if the fields number is
        %greater than the number of corners. It uses the fields with max
        %cornerscore_ea minus the extra fields cornerscorep (penalty score)
        ncorners = size(env_coor,1)-1;
        if penalty
            if mask
                %first, to penalize the cornerscore inside the mask, if
                %the number of fields in the mask greater than the number
                %of corners. 
                if numel(cornerscore_ea) > ncorners
                    [maxcs,maxcs_idx] = maxk(cornerscore_ea, ncorners);
                    idx = ~ismember([1:numel(cornerscore_ea)],maxcs_idx);
                    cornerscorep = abs(cornerscore_ea(idx)-1); %to make penalty linear
                    cornerscore = (sum(maxcs) - sum(cornerscorep))/ncorners;
                end
                %second, get the ratemap outside the mask. If there are
                %fields outside the mask, panelize the cornerscore with the
                %score for the outside fields. 
                BW = imregionalmax(rtmap);
                BW(msk) = 0;
                [rowm,colm] = find(BW == 1);
                if ~isempty(rowm)
                    field_coorm = [colm,rowm];
                    Dm = pdist2(env_coor,field_coorm);
                    cornerscorep2 = (Dm(end,:) - min(Dm(1:end-1,:)))./(Dm(end,:) + min(Dm(1:end-1,:)));
                    cornerscorep2 = abs(cornerscorep2-1); 
                    if ~exist('cornerscore','var')
                        cornerscore = (sum(cornerscore_ea) - sum(cornerscorep2))/ncorners;
                    else
                        cornerscore = cornerscore - sum(cornerscorep2)/ncorners;
                    end
                else
                    if ~exist('cornerscore','var')
                        cornerscore = sum(cornerscore_ea)/ncorners;
                    end
                end
            else
                if numel(cornerscore_ea) > ncorners
                    [maxcs,maxcs_idx] = maxk(cornerscore_ea, ncorners);
                    idx = ~ismember([1:numel(cornerscore_ea)],maxcs_idx);
                    cornerscorep = abs(cornerscore_ea(idx)-1); %to make penalty linear
                    cornerscore = (sum(maxcs) - sum(cornerscorep))/ncorners;
                else
                    cornerscore = sum(cornerscore_ea)/ncorners;
                end
            end
        else
            cornerscore = sum(cornerscore_ea)/ncorners;
        end
        %allocate the final results
        cmetrics.field_coor{ii} = field_coor;
        cmetrics.field2corner_dist{ii} = D;
        cmetrics.cornerscore_ea{ii} = cornerscore_ea;
        if exist('cornerscorep','var')
            cmetrics.cornerscorep{ii} = cornerscorep;
        end
        cmetrics.cornerscore(ii) = cornerscore;
        %clear variables
        if exist('cornerscorep','var')
            clear cornerscorep
        end
        if exist('cornerscore','var')
            clear cornerscore
        end
    end
    
end

end

