function [num_motif,motif_idx,num_iso_cluster] = extract_motif(S,p,min_burst)
if nargin<2, p=mean(S(:)); end
if nargin<3, min_burst = sqrt(size(S,1)); end

[idx,netsim,dpsim,expref]=apcluster(S,p,'maxits',1000,'convits',100,'dampfact',0.7,'nonoise');

uidx = unique(idx);

motif_idx = cell(1,1);
num_motif = 0;
num_iso_cluster = 0;
for i = 1:numel(uidx)
    ix = find(idx==uidx(i));
    
    if numel(ix)>min_burst
        num_motif = num_motif + 1;
        motif_idx{num_motif} = ix;
    end
    
    if numel(ix)<2
        num_iso_cluster = num_iso_cluster+1;
    end
end