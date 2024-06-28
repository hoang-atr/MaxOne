function Cmat = construct_connectivity_matrix(spike_time, temporal_window)
if nargin<2, temporal_window=1; end

Nelectrode = numel(spike_time);

Cmat = zeros(Nelectrode,Nelectrode);
for i = 1:Nelectrode-1    
    for j = (i+1):Nelectrode       
            spk1 = spike_time{i};
            spk2 = spike_time{j};
            n1 = numel(spk1);
            n2 = numel(spk2);

            if (n1==0) || (n2==0), continue, end

            for ii = 1:n1
                d = abs(spk2 - spk1(ii));
                ix = d<temporal_window;
                Cmat(i,j) = Cmat(i,j) + sum(ix);
            end   
            Cmat(i,j) = Cmat(i,j) ./ sqrt(n1*n2);       
    end
end