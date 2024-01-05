function outputMaps = gaussianfilter2D(inputMaps, halfNarrow, narrowStdev)
%This is a 2D Gaussian filter to smooth 2D ratemap smoothing function
%It's original form was obtain from Dr. Doug Nitz at UCSD
%Modified by Yanjun Sun, 12/28/22

if ~exist('narrowStdev', 'var') || isempty(narrowStdev)
    narrowStdev = 2;
end
if ~exist('halfNarrow', 'var') || isempty(halfNarrow)
    halfNarrow = 3; %when binsize==2
end

[xGridVals, yGridVals]=meshgrid(-halfNarrow:1:halfNarrow);
narrowGaussian = exp(-0.5 .* (xGridVals.^2+yGridVals.^2)/narrowStdev^2)/(narrowStdev*(sqrt(2*pi)));
narrowGaussianNormed=narrowGaussian./sum(sum(narrowGaussian));

for i = 1:size(inputMaps,3)
    outputMaps(:,:,i) = nanconv(inputMaps(:,:,i,1),narrowGaussianNormed, 'nanout');
end

end


