function neuron = neuronPlot(neuron,frames_pass)

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
    plot(neuron.trace(i,frames_pass),'k')
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
    plot(neuron.S(i,frames_pass),'r')
    set(gca,'Xtick',[]);set(gca,'FontSize',8)
    axis tight
    if i == 1,title('Filtered data');end
end

% correct the data
neuron.C = neuron.C(:,frames_pass);
neuron.C_raw = neuron.C_raw(:,frames_pass);
neuron.C_df = neuron.C_df(:,frames_pass);
neuron.S = neuron.S(:,frames_pass);
neuron.trace = neuron.trace(:,frames_pass);
neuron.time = neuron.time(frames_pass);
neuron.num2read = length(frames_pass);
