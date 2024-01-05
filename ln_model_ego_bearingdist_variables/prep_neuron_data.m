function [neuronIndiv, behavIndiv] = prep_neuron_data(neuronIndiv_ori,...
    behavIndiv_ori, c_vector)
% to combine data from different sessions
% c_vector is a vector or matrix that is same size as orginal data but
% specified the to-be combined sessions with same values, see below

if ~exist('c_vector','var')||isempty(c_vector)
    c_vector = [1,1,2,2;3,3,4,4];
end
    matcol = max(max(c_vector))/size(c_vector,1);%col of cell array after combination
    neuronIndiv = cell(1,max(max(c_vector)));
    behavIndiv = cell(1,max(max(c_vector)));
for ii = 1:max(max(c_vector))
    [row,col] = find(c_vector == ii);
    [N,B] = combine_sessions(neuronIndiv_ori,behavIndiv_ori,row,col);
    neuronIndiv{1,ii} = N;
    behavIndiv{1,ii} = B;
end
    %reshape if orignal data has two rows, hard coded
    if size(neuronIndiv,2) > matcol
        neuronIndiv = reshape(neuronIndiv,matcol,2)';
        behavIndiv = reshape(behavIndiv,matcol,2)';
    end
    
end