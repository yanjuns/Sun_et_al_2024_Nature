function de_corners = compare_decodeError_at_corners(decodeResults,decodeResultsKO)
%This function aims to find the decoding errors at each phyical corner.
%the last row of each elements is the decoding error in the center of the
%arena. 
%yanjuns@stanford.edu, 12/27/22

load('env_geometry.mat','env_coori')
envc = env_coori;
de_corners = cell(1,length(envc)); %pkfr at each corner
for n = 1:length(envc)
    envc_ea = envc{1,n};
    de2 = [];de2KO = [];
    for jj = 1:size(decodeResults,1)
        de2 = cat(3,de2, decodeResults{jj,n}.decodeError2_2step);
        de2KO = cat(3,de2KO, decodeResultsKO{jj,n}.decodeError2_2step);
    end
    de2 = nanmean(de2,3);
    de2KO = nanmean(de2KO,3);
    de = de2KO./de2;
%     de = gaussianfilter2D(de, 3, 1.5);

    v = ones(8);
    de_corners_ea = [];
    for ii = 1:size(envc_ea,1)
        M = zeros(size(de,1),size(de,2));
        M(envc_ea(ii,2),envc_ea(ii,1)) = 1;
        corner_regidx = logical(conv2(M,v,'same'));
        corner_region = de(corner_regidx);
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

