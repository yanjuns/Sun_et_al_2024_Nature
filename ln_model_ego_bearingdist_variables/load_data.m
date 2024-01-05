%% load_data
switch p.dataType
    case{'whole','Whole'}
        load('neuronIndividualsf.mat','neuronIndividualsf')
        load('behavIndividualsf.mat','behavIndividualsf')
        load('thresh.mat','thresh')
        load('env_geometry.mat','S')
        if p.c
            [neuronIndiv, behavIndiv] = prep_neuron_data(neuronIndividualsf,...
                behavIndividualsf, p.c_vector);
        else
            neuronIndiv = neuronIndividualsf;
            behavIndiv = behavIndividualsf;
        end
        
    case{'subset','Subset','Separated','separated','LR','divided','Divided'}
        load('NBindivLR.mat','neuronIndivLR','behavIndivLR')
        load('thresh.mat','thresh')
        load('env_geometry.mat','S')
        if p.c
            [neuronIndiv, behavIndiv] = prep_neuron_data(neuronIndivLR,...
                behavIndivLR, p.c_vector);
        else
            neuronIndiv = neuronIndivLR;
            behavIndiv = behavIndivLR;
        end
end