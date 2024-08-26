

function [update_cluster_disks_pos] = BoundaryCondition_clusters(update_cluster_disks_pos, XX, YY);

cluster_x = update_cluster_disks_pos(:,2);
cluster_y = update_cluster_disks_pos(:,3);

mean_xc = mean(cluster_x);
mean_yc = mean(cluster_y);

% shifting positions of all disks in a cluster based on position of center of mass of the cluster

if mean_xc < 0
    cluster_x = cluster_x + XX;
end

if mean_xc > XX
    cluster_x = cluster_x - XX;
end

if mean_yc < 0
    cluster_y = cluster_y + YY;
end

if mean_yc > YY
    cluster_y = cluster_y - YY;
end

update_cluster_disks_pos(:,2) = cluster_x;
update_cluster_disks_pos(:,3) = cluster_y;

end