
function clusters_ids = Clusters_finding_new(disks_stat, bonded_pairs_sorted_unique)

num_particles = size(disks_stat, 1);
bonds = bonded_pairs_sorted_unique;
adj_list = cell(num_particles, 1);

for i = 1:size(bonds, 1)
    adj_list{bonds(i, 1)} = [adj_list{bonds(i, 1)}, bonds(i, 2)];
    adj_list{bonds(i, 2)} = [adj_list{bonds(i, 2)}, bonds(i, 1)];
end

visited = false(num_particles, 1);
clusters_ids = {};

for i = 1:num_particles
    if ~visited(i)
       cluster = [];
       cluster = dfs(i, cluster);
       clusters_ids{end+1} = cluster;
    end
end

function cluster = dfs(p, cluster)
 visited(p) = true;
 cluster(end+1) = p; % Add particle to current cluster
 for neighbor = adj_list{p}
     if ~visited(neighbor)
        cluster = dfs(neighbor, cluster); % Recursively visit neighbors
     end
 end
end

clusters_ids = clusters_ids';
clusters_ids_sorted = {};
for i = 1:length(clusters_ids)
    aa = clusters_ids{i,1};
    aa = aa';
    if length(aa) ~= 1
       clusters_ids_sorted = [clusters_ids_sorted; aa];
    end
end

clusters_ids = clusters_ids_sorted;


end


