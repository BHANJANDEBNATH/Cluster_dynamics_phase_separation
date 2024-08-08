
% overlap check for clusters

function [condition_cluster] = overlap_check_clusters(ids, update_cluster_disks_pos, disks_stat, dia_disk, XX, YY)

condition_cluster = [];

for i = 1:length(ids)
    posx = update_cluster_disks_pos(i,2);
    posy = update_cluster_disks_pos(i,3);
    
    % for overlap check among all disks, it is mandatory that all disks to be
    % inside the simulation box
    % it may happen that though the center of mass of cluster is inside simulation
    % box, but some disks corrsponding to that cluster are outside the simulation box 
    % because we are treating each cluster (say,) as one large entity
    % note the center of mass of cluster is always inside the simulation box (see overlap_check_clusters.m)
    % that's why, shifting positions of those disks in the cluster outside
    % simulation domain is mandatory to check overlap with other disks
    
    if posx <= 0
        posx = posx + XX;
    end
    if posx > XX
        posx = posx - XX;
    end
    if posy <= 0
        posy = posy + YY;
    end
    if posy > YY
        posy = posy - YY;
    end

    for j = 1:length(disks_stat)
        if ids(i,1) ~= disks_stat(j,1)
            dis = ((posx - disks_stat(j,2))^2 + (posy - disks_stat(j,3))^2)^0.5;
            if dis > dia_disk
                condition_cluster = [condition_cluster; 0];
            else
                condition_cluster = [condition_cluster; 1];
            end
        end
    end
end


end