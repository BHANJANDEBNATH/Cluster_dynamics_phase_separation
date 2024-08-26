
% Brownian dynamics of clusters
% moving clusters based on center of mass

function[disks_stat] = BD_Clusters_pos_update(ids,cluster_disks_pos_old,disks_stat,D_free_disk,dia_disk,tau,XX,YY)


% calculate radius of gyration of cluster
        xcom = 0;
        ycom = 0;
        for mm = 1:length(ids)
              xcom = xcom + cluster_disks_pos_old(mm,2);
              ycom = ycom + cluster_disks_pos_old(mm,3);
        end
        xcom = xcom/length(ids);
        ycom = ycom/length(ids);
        xr = 0;
        yr = 0;
        for mm = 1:length(ids)
            xr = xr + (cluster_disks_pos_old(mm,2) - xcom)^2;
            yr = yr + (cluster_disks_pos_old(mm,3) - ycom)^2;
        end
        rad_gyr = (xr + yr)^0.5/length(ids);

        % translational and rotational difusion coefficients of cluster
        tr_diff_cluster = D_free_disk * (dia_disk/2)/rad_gyr;
        rot_diff_cluster = (6/8) * D_free_disk * (dia_disk/2) * (1/rad_gyr)^3;

        % translate and rotate cluster with respect to center of mass
        dis = sqrt(2 * 1 * tr_diff_cluster * tau); 
        ds = dis * randn(1,2);                          
        theta_rad = sqrt(2 * 1 * rot_diff_cluster * tau) * randn(1,1);  % in radian
        [update_cluster_disks_pos] = translate_rotate_cluster(cluster_disks_pos_old, ds, theta_rad);

        % periodic BC for clusters
        [update_cluster_disks_pos] = BoundaryCondition_clusters(update_cluster_disks_pos, XX, YY);

        % check overlap for disks of clusters with other clsters and disks
        [condition_cluster] = overlap_check_clusters(ids, update_cluster_disks_pos, disks_stat, dia_disk, XX, YY);
        
        if mean(condition_cluster) == 0
           [r,c] = size(update_cluster_disks_pos);
           for kk = 1:r
               iii = update_cluster_disks_pos(kk,1);
               disks_stat(iii,2:3) = update_cluster_disks_pos(kk,2:3);
           end
        end



end



