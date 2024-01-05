function de_corners = find_decodeError_at_corners(decodeResults)
%This function aims to find the decoding errors at each phyical corner.
%the last row of each elements is the decoding error in the center of the
%arena. 
%yanjuns@stanford.edu, 12/27/22

load('env_geometry.mat','env_coori','S')
envc = env_coori;
de_corners = cell(1,length(envc)); %pkfr at each corner
for n = 1:length(envc)
    envc_ea = envc{1,n};
    de2 = [];
    for jj = 1:size(decodeResults,1)
        de2 = cat(3,de2, decodeResults{jj,n}.decodeError2_2step);
    end
    de2 = nanmean(de2,3);
    v = ones(12);
    de_corners_ea = [];
    for ii = 1:size(envc_ea,1)
        M = zeros(size(de2,1),size(de2,2));
        M(envc_ea(ii,2),envc_ea(ii,1)) = 1;
        corner_regidx = logical(conv2(M,v,'same'));
        corner_region = de2(corner_regidx);
        de_corners_ea = [de_corners_ea;nanmean(nanmean(corner_region))];
%         %to plot each extracted region
%         row = sum(any(corner_regidx > 0, 2));
%         col = sum(any(corner_regidx > 0, 1));
%         corner_region = reshape(corner_region, [row,col]);
%         figure
%         imagesc(corner_region)
%         axis image
    end
    de_corners{n} = de_corners_ea;
end

end

