function [neuronL, neuronR, behavL, behavR, autocenter] = split_neuron_behav_LR(neuron, behav, centerdetect, edgeCorrect)
% This function aims to split CPP neuron and behav data under baseline or
% test condition into left and right sub-datasets and facilitate downstream
% data analysis.
% Yanjun Sun, Stanford University, 9/10/2019
if ~exist('edgeCorrect','var') || isempty(edgeCorrect)
    edgeCorrect = true;
end
if ~exist('centerdetect','var') || isempty(centerdetect)
    centerdetect = 'auto';
end

neuronL = struct; neuronR = struct; 
behavL = struct; behavR = struct;
%% detect the border between two CPP chambers
ymax = max(behav.position(:,2));
ymin = min(behav.position(:,2));
xidx = find(behav.position(:,2)>0 & behav.position(:,2)< (ymax+ymin)/10*3 |...
    behav.position(:,2)>ymax/10*7 & behav.position(:,2)< ymax);
xban = behav.position(xidx,:);
switch centerdetect
    case{'fixed', 'hard', 'physical'}
        %simple physical center
        autocenter = (max(behav.position(:,1))+min(behav.position(:,1)))/2;
    case{'auto'}
        %auto center detection by using chamber patterns
        binedge = ceil(min(xban(:,1))):0.1:floor(max(xban(:,1)));
        [N,~] = histcounts(xban(:,1),binedge);
        idx = find(N==0);
        autocenter = median(binedge(idx));
    case{'mannual','enter','mannually enter','mannual enter'}
        figure;
        plot(behav.position(:,1), behav.position(:,2));
        axis image
%         answer = inputdlg('Enter the CPP center value after examine the figure:',...
%              'CPP box center', [1 60])
%         autocenter = str2double(answer);
        prompt = 'Enter the CPP center value after examine the figure:';
        autocenter = input(prompt)
        close(gcf);

end
figure;
plot(behav.position(:,1), behav.position(:,2));
axis image
hold on;
y = linspace(ceil(min(behav.position(:,2))),floor(max(behav.position(:,2))));
x = repelem(autocenter, 100);
plot(x, y, 'r');
%% results output
if edgeCorrect
    % edge correction parameters
    xbanright = xban(xban(:,1) > autocenter);
    edgediff = min(behav.position(:,1)) - (min(xbanright(:,1))-autocenter);
    % behavL and behavR
    neuron.pos = interp1(behav.time,behav.position,neuron.time);
    idxbposL = find(behav.position(:,1) < autocenter);
    idxbposR = find(behav.position(:,1) > autocenter);
    behavL.position = behav.position(idxbposL,:);
    behavL.positionblue = behav.positionblue(idxbposL,:);
    behavL.time = behav.time(idxbposL,:);
    behavL.speed = behav.speed(:,idxbposL);
    behavL.hd = behav.hd(idxbposL,:);
    behavR.position = behav.position(idxbposR,:);
    behavR.position(:,1) = behavR.position(:,1) - autocenter + edgediff; %min(behavL.position(:,1)); % edge correction
    behavR.positionblue = behav.positionblue(idxbposR,:);
    behavR.positionblue(:,1) = behavR.positionblue(:,1) - autocenter + edgediff;
    behavR.time = behav.time(idxbposR,:);
    behavR.speed = behav.speed(:,idxbposR);
    behavR.hd = behav.hd(idxbposR,:);
    %neuronL and neuronR
    idxnposL = find(neuron.pos(:,1) < autocenter);
    idxnposR = find(neuron.pos(:,1) > autocenter);
    neuronL.pos = neuron.pos(idxnposL,:);
    neuronL.time = neuron.time(idxnposL,:);
    neuronL.S = neuron.S(:,idxnposL);
    neuronL.trace = neuron.trace(:,idxnposL);
    neuronR.pos = neuron.pos(idxnposR,:);
    neuronR.pos(:,1) = neuronR.pos(:,1)- autocenter + edgediff; %min(behavL.position(:,1)); % edge correction
    neuronR.time = neuron.time(idxnposR,:);
    neuronR.S = neuron.S(:,idxnposR);
    neuronR.trace = neuron.trace(:,idxnposR);
else
    % behavL and behavR
    neuron.pos = interp1(behav.time,behav.position,neuron.time);
    idxbposL = find(behav.position(:,1) < autocenter);
    idxbposR = find(behav.position(:,1) > autocenter);
    behavL.position = behav.position(idxbposL,:);
    behavL.positionblue = behav.positionblue(idxbposL,:);
    behavL.time = behav.time(idxbposL,:);
    behavL.speed = behav.speed(:,idxbposL);
    behavL.hd = behav.hd(idxbposL,:);
    behavR.position = behav.position(idxbposR,:);
    behavR.positionblue = behav.positionblue(idxbposR,:);
    behavR.time = behav.time(idxbposR,:);
    behavR.speed = behav.speed(:,idxbposR);
    behavR.hd = behav.hd(idxbposR,:);
    %neuronL and neuronR
    idxnposL = find(neuron.pos(:,1) < autocenter);
    idxnposR = find(neuron.pos(:,1) > autocenter);
    neuronL.pos = neuron.pos(idxnposL,:);
    neuronL.time = neuron.time(idxnposL,:);
    neuronL.S = neuron.S(:,idxnposL);
    neuronL.trace = neuron.trace(:,idxnposL);
    neuronR.pos = neuron.pos(idxnposR,:);
    neuronR.time = neuron.time(idxnposR,:);
    neuronR.S = neuron.S(:,idxnposR);
    neuronR.trace = neuron.trace(:,idxnposR);
end
end