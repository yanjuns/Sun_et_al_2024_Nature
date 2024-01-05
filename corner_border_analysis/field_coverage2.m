function coverage = field_coverage2(fields)
    coverage.perimeter.W=0;coverage.perimeter.N=0;coverage.perimeter.E=0;coverage.perimeter.S=0;coverage.area.total=0;
    coverage.area.inside_area=0;coverage.area.W=0;coverage.area.E=0;coverage.area.S=0;
    coverage.area.N=0; coverage.weighted_distance=0;
if(length(fields)>0)
    [ly,lx]=size(fields{1});
    l=round(min(lx,ly)/4); 
    coverage.perimeter.max_one_field=0;
%     total_perimeter=0;
    for i=1:length(fields)
        
        aux_map=fields{i}(:,1:8);
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=covered/norm; end 
            
        aux_map=fields{i}(:,end:-1:end+1-8);       
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=covered/norm; end 
        
        aux_map=fields{i}(1:8,:)'; 
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=covered/norm; end         

        aux_map=fields{i}(end:-1:end+1-8,:)'; 
        [covered,norm]=wall_field(aux_map);
        if(covered/norm>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=covered/norm; end                 
        
        aux=sum(fields{i}(1:end,1)>0);
        %if(aux/ly>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=aux/ly; end 
        coverage.perimeter.W=coverage.perimeter.W+aux;
        aux=sum(fields{i}(1:end,end)>0);
        %if(aux/ly>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=aux/ly;end
        coverage.perimeter.E=coverage.perimeter.E+aux;
        aux=sum(fields{i}(1,1:end)>0);
%        if(aux/lx>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=aux/lx; end
        coverage.perimeter.S=coverage.perimeter.S+aux;
        aux=sum(fields{i}(end,1:end)>0);
%        if(aux/lx>coverage.perimeter.max_one_field) coverage.perimeter.max_one_field=aux/lx; end
        coverage.perimeter.N=coverage.perimeter.N+aux;
        coverage.area.total=coverage.area.total+sum(fields{i}(:)>0);
        coverage.area.inside_area=coverage.area.inside_area+sum(sum(fields{i}(l:ly+1-l,l:lx+1-l)>0));
        coverage.area.W=coverage.area.W+sum(sum(fields{i}(1:end,1:l)>0));
        coverage.area.E=coverage.area.E+sum(sum(fields{i}(1:end,lx+1-l:lx)>0));
        coverage.area.S=coverage.area.S+sum(sum(fields{i}(1:l,1:end)>0));
        coverage.area.N=coverage.area.N+sum(sum(fields{i}(ly+1-l:end,1:end)>0)); 
        
        map_rot=fields{i};
        cov_aux=search_peak_distance(map_rot);
        map_rot=flipud(map_rot);
        cov_aux=[cov_aux search_peak_distance(map_rot)];
        map_rot=fliplr(map_rot);
        cov_aux=[cov_aux search_peak_distance(map_rot)];
        map_rot=flipud(map_rot);
        cov_aux=[cov_aux search_peak_distance(map_rot)];  
        coverage.distance_peak_to_wall{i}=nanmean(cov_aux);
    end
    coverage.area.outside_area=coverage.area.total-coverage.area.inside_area;
    coverage.area.inside_area=coverage.area.inside_area/((ly+2-2*l)*(lx+2-2*l));
    coverage.area.outside_area=coverage.area.outside_area/(ly*lx-(ly+2-2*l)*(lx+2-2*l));
    coverage.area.W=coverage.area.W/ly/l;
    coverage.area.E=coverage.area.E/ly/l;
    coverage.area.S=coverage.area.S/lx/l;
    coverage.area.N=coverage.area.N/lx/l;
    
    [coverage.area.max,walli] = max([coverage.area.N,coverage.area.S,coverage.area.E,coverage.area.W]);
    w = {'N','S','E','W'};
    coverage.area.wall = w(walli);   
    coverage.perimeter.total=coverage.perimeter.W+coverage.perimeter.E+coverage.perimeter.N+coverage.perimeter.S;
    coverage.perimeter.total=coverage.perimeter.total/(2*lx+2*ly);
    coverage.perimeter.W=coverage.perimeter.W/ly;
    coverage.perimeter.E=coverage.perimeter.E/ly;    
    coverage.perimeter.S=coverage.perimeter.S/lx;    
    coverage.perimeter.N=coverage.perimeter.N/lx;  
    [coverage.perimeter.max,walli] = max([coverage.perimeter.N,coverage.perimeter.S,coverage.perimeter.E,coverage.perimeter.W]);
    coverage.perimeter.wall = w(walli);
    coverage.area.total=coverage.area.total/lx/ly;
end
function [covered,norm]=wall_field(map)
[ly,~]=size(map);
for j=1:ly
    a=find(isfinite(map(j,:)),1,'first');
    if length(a)>0
        aux(j)=map(j,a);
    else
        aux(j)=nan;
    end
end
norm=sum(isfinite(aux));
covered=nansum(aux>0);
function distance_vector=search_peak_distance(map)
distance_vector=[];
aux=find(map(1,:)>0);
aux_max=min([1:size(map,2); size(map,2):-1:1]);
halfw = floor(min(size(map,1),size(map,2))/2);
for i=1:length(aux)
    j=aux(i);
    maxL = aux_max(j);
    if maxL > halfw
        maxL = halfw;
    end
    peaks=find_peaks_smooth(map(1:maxL,aux(i))',5,0);
    if length(peaks)>0
        distance_vector(i)=peaks(1);
    else
        distance_vector(i)=nan;
    end
end