function [A,As,pA,pAs] = anat_density_map(C,celltype,binsize)
%This function aims to generate the anatomical density map for defined
%functional cell types. 
%yanjuns@stanford.edu, 3/31/2023

%calculate density for all neurons
[A,As,xAxis,yAxis] = anat_density(C, [], [], 1, binsize);

%place cell density (normalized by all neurons)
pC = C(celltype,:);
[pN,pNs,~,~] = anat_density(pC, xAxis, yAxis);
pA = pN./A;
pA(isnan(pA)) = 0;
pAs = imgaussfilt(pA,1);

end

