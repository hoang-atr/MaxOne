function [percentile, cor] = find_optimal_electrode(spike_time, duration, binsize)
Nelectrode = numel(spike_time);

count_spike = cellfun(@numel, spike_time);

edge=0:binsize:duration;
Nt = numel(edge)-1;

tsr_base = zeros(1, Nt);
for n = 1:Nelectrode
    count = histcounts(spike_time{n}, edge);                            
    count(count>1)=1;
    tsr_base = tsr_base + count;
end

percentile = 0.1:0.1:0.9;
Np = numel(percentile);

cor = zeros(1,Np);
for p = 1:numel(percentile)
    idx = count_spike>quantile(count_spike,percentile(p));
    spike_perc = spike_time(idx);
    Nelectrode_perc = numel(spike_perc);
    
    tsr_perc = zeros(1, Nt);
    for n = 1:Nelectrode_perc
        count = histcounts(spike_perc{n}, edge);                            
        count(count>1)=1;
        tsr_perc = tsr_perc + count;
    end
    
    idxNonZero = tsr_base>0 | tsr_perc>0;
    tmp = corrcoef(tsr_base(idxNonZero), tsr_perc(idxNonZero));
    cor(p) = tmp(1,2);        
end


