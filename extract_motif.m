function motif_idx = extract_motif(S,p,min_burst)
%% input 
% S: - similarity matrix; 
% p:  preference parameter of Affinity Propagation clustering
% min_burst: minimum burst size of clusters to be considered as motifs
%% output
% motif_idx: index of bursts that belong to the repeating motifs
% Note: size of <motif_idx> is the number of repeating motifs
Nburst = size(S,1);

if nargin<2, p=median(S(:)); end
if nargin<3, min_burst = sqrt(Nburst); end

%% Affinity Propagation clustering
idx=apcluster(S,p,'maxits',1000,'convits',100,'dampfact',0.7,'nonoise');

%% check the cluster size
uidx = unique(idx);
cluster_idx = cell(1,1);
num_cluster = 0;
for i = 1:numel(uidx)
    ix = find(idx==uidx(i));
    
    if numel(ix)>min_burst
        num_cluster = num_cluster + 1;
        cluster_idx{num_cluster} = ix;
    end    
end

%% check the separation between clusters
% similarity within-cluster
Swi = cell(num_cluster,1);
for i = 1:num_cluster
    idx = cluster_idx{i};
    Sm = S(idx,idx);
    ix = triu(Sm,1)>eps; 
    
    Swi{i} = Sm(ix);
end

% check the separation between within- and across-cluster
motif = ones(num_cluster,1);
for i = 1:num_cluster
    idx = cluster_idx{i};
    idx_k = setdiff(1:Nburst,idx);
    Sac = S(idx,idx_k); % similarity across-cluster
    [~,p] = ttest2(Sac(:),Swi{i},'tail','right');
    if p < 0.01
        motif(i) = 0;        
    end
end

motif_idx=cluster_idx(motif>0);
