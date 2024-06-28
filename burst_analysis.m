function [tau_act,tau_st,S,A,first_spike_in_burst,burst_first_t,burst_end_t] = burst_analysis(spike_time, binsize, sd_tsr, time_threshold, min_burst_duration)
if nargin<2, binsize = 10; end
if nargin<3, sd_tsr = 0.2; end
if nargin<4, time_threshold = 50; end
if nargin<5, min_burst_duration = 50; end

Nelectrode = numel(spike_time);

duration = 0;
for n = 1:Nelectrode
    tmp = max(spike_time{n});
    if duration < tmp
        duration = tmp;
    end
end
duration = duration + 1;

edge=0:binsize:duration;
Nt = numel(edge)-1;

tsr = zeros(1, Nt);
for n = 1:Nelectrode
    count = histcounts(spike_time{n}, edge); 
    count(count>1) = 1;
    tsr = tsr + count;
end

burst_theshold = sd_tsr*std(tsr);
burst_bin = find(tsr>burst_theshold);

d = diff([0,diff(burst_bin)==1,0]);
burst_first_t = edge(burst_bin(d>0))-binsize;
burst_end_t = edge(burst_bin(d<0))+binsize;
Nburst = numel(burst_first_t);

first_spike_in_burst = zeros(Nelectrode, Nburst);
last_spike_in_burst = zeros(Nelectrode, Nburst);
burst_duration = zeros(Nburst,1);
all_spikes_in_burst = cell(Nelectrode, Nburst);

for b = 1:Nburst
    for n = 1:Nelectrode
        spk = spike_time{n};
        ix = find(spk>=burst_first_t(b) & spk<=burst_end_t(b));
        
        all_spikes_in_burst{n,b} = spk(ix);
        
        if numel(ix)==0
            first_spike_in_burst(n,b) = nan;
            last_spike_in_burst(n,b) = nan;
        else
            first_spike_in_burst(n,b) = spk(ix(1));
            last_spike_in_burst(n,b) = spk(ix(end));
        end
    end
    burst_duration(b) = max(last_spike_in_burst(:,b))-min(first_spike_in_burst(:,b));
end

idx = burst_duration<min_burst_duration;
burst_duration(idx) = []; 
burst_first_t(idx) = [];
burst_end_t(idx) = [];
first_spike_in_burst(:,idx)=[];
last_spike_in_burst(:,idx)=[];
Nburst = numel(burst_duration);

tau_act = max(first_spike_in_burst)-min(first_spike_in_burst);
tau_st = mean(last_spike_in_burst-first_spike_in_burst,'omitnan');

A = cell(Nburst,1);
for b = 1:Nburst
    fs = first_spike_in_burst(:,b);
    a = zeros(Nelectrode,Nelectrode);
    for i = 1:Nelectrode
        for j = 1:Nelectrode
            a(i,j) = fs(i)-fs(j);
        end
    end
    A{b} = a;    
end

S = zeros(Nburst,Nburst);
for p = 1:Nburst
    for q = 1:Nburst
        da = abs(A{p} - A{q});
        not_diag_idx = eye(Nelectrode)<eps;
        ix = da(not_diag_idx)<=time_threshold & ~isnan(da(not_diag_idx));
        S(p,q) = sum(ix(:))/(Nelectrode*(Nelectrode-1));
    end
end
