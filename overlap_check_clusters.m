
% overlap check for clusters

function [condition_cluster] = overlap_check_clusters(ids,update_cluster_disks_pos,disks_stat,dia_disk,neighborList)

condition_cluster = [];

for i = 1:length(ids)
    imp = ids(i,1);
    posx = update_cluster_disks_pos(i,2);
    posy = update_cluster_disks_pos(i,3);
    
    % sc = disks_stat(:,1:3);   % for Brute-Force check
    
    neighbor_ids = neighborList{imp,1};
    neighbor_ids = neighbor_ids';
    sc = [];
    for m = 1:length(neighbor_ids)
        neigh = neighbor_ids(m,1);
        sc = [sc; disks_stat(neigh,1:3)];
    end


    
    [r c] = size(sc);
    for j = 1:r
        if imp ~= sc(j,1)
            xx2 = sc(j,2);
            yy2 = sc(j,3);
            % distance
            dis = ((posx - xx2)^2 + (posy - yy2)^2)^0.5;
            if dis > dia_disk
                condition_cluster = [condition_cluster; 0];
            else
                condition_cluster = [condition_cluster; 1];
            end
        end
    end

end


end