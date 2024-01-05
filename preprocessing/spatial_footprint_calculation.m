function [footprint]=spatial_footprint_calculation(neuron,varargin)

d1=neuron.options.d1;
d2=neuron.options.d2;
summ=zeros(d1*d2,1);

if length(varargin)>0
    thresh=varargin{1};
else
    thresh=0;
end

for i=1:size(neuron.A,2)
    Ai=neuron.A(:,i);
    Ai=Ai*1/max(Ai(:));
    ind=find(Ai>thresh);%the pixels at which a neuron shows up
    summ(ind)=Ai(ind);
end

footprint=reshape(summ,d1,d2);
