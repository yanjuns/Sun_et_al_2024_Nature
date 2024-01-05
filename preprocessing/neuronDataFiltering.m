function [neuron,frame_pass_neuron] = neuronDataFiltering(neuron,frame_pass,behav)
%% calculate the corresponding frames in neuron data based on the filtered frames in behav data
% find the nearest points based on the time
frame_filter_behav = setdiff(1:length(behav.time),frame_pass);
frame_filter_neuron = knnsearch(neuron.time,behav.time(frame_filter_behav));
frame_pass_neuron = setdiff(1:size(neuron.C,2),frame_filter_neuron);

%% plotting
figure
ha = tight_subplot(5,2,[.01 .04],[.02 .03],[.2 .2]);
for i = 1:5
    axes(ha(2*i-1))
    % subplot(2,2,1)
    plot(neuron.trace(i,:),'k')
    hold on
    %     plot(neuron.S(i,:),'r')
    set(gca,'Xtick',[]);set(gca,'FontSize',8)
    ylabel(['Neuron ',num2str(i)],'FontSize',10)
    axis tight
    if i == 1,title('Original data');end
    
    % subplot(2,2,2)
    axes(ha(2*i))
    plot(neuron.trace(i,frame_pass_neuron),'k')
    hold on
    %     plot(neuron.S(i,frames_pass),'r')
    set(gca,'Xtick',[]);set(gca,'FontSize',8)
    axis tight
    if i == 1,title('Filtered data');end
end

figure
ha = tight_subplot(5,2,[.01 .04],[.02 .03],[.2 .2]);
for i = 1:5
    axes(ha(2*i-1))
    % subplot(2,2,1)
    plot(neuron.S(i,:),'r')
    set(gca,'Xtick',[]);set(gca,'FontSize',8)
    ylabel(['Neuron ',num2str(i)],'FontSize',10)
    axis tight
    if i == 1,title('Original data');end
    
    % subplot(2,2,2)
    axes(ha(2*i))
    plot(neuron.S(i,frame_pass_neuron),'r')
    set(gca,'Xtick',[]);set(gca,'FontSize',8)
    axis tight
    if i == 1,title('Filtered data');end
end

% correct the data
neuron.C = neuron.C(:,frame_pass_neuron);
neuron.C_raw = neuron.C_raw(:,frame_pass_neuron);
% neuron.C_df = neuron.C_df(:,frame_pass_neuron);
neuron.S = neuron.S(:,frame_pass_neuron);
neuron.trace = neuron.trace(:,frame_pass_neuron);
neuron.trace_raw = neuron.trace_raw(:,frame_pass_neuron);
neuron.time = neuron.time(frame_pass_neuron);
neuron.num2read = length(frame_pass_neuron);
