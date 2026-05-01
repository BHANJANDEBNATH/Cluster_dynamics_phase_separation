
% Brownian dynamics of clusters
% moving clusters based on center of mass

function[disks_stat] = BD_Clusters_pos_update(ids,cluster_disks_pos_old,disks_stat,D_free_disk,dia_disk,tau,XX,YY,neighborList)

        positions = cluster_disks_pos_old(:,2:3);
        [rad_gyr,unwrapped_positions] = radius_gyration_pbc(ids,positions,XX,YY, dia_disk);

        % translational and rotational difusion coefficients of cluster
        tr_diff_cluster = D_free_disk * (dia_disk/2)/rad_gyr;

        rot_diff_cluster = (1/2) * tr_diff_cluster  * (1/rad_gyr)^2;
        
        % translate and rotate cluster
        dis_c = sqrt(2 * 1 * tr_diff_cluster * tau); 
        ds = dis_c * randn(1,2);                          
        theta_rad = sqrt(2 * 1 * rot_diff_cluster * tau) * randn(1,1);      % in radian
        theta_rad = mod(theta_rad, 2 * pi);                                 % normalize angle to [0 2*pi]
      
        translate = [ds(1,1) ds(1,2)];
        rotate = theta_rad;
        pbc_xmin = 0; 
        pbc_xmax = XX;
        pbc_ymin = 0; 
        pbc_ymax = YY;
       
        % translate and rotate 
        [update_coord] = translate_rotate_cluster(ids,unwrapped_positions,translate,rotate,pbc_xmin,pbc_xmax,pbc_ymin,pbc_ymax);

        update_cluster_disks_pos = cluster_disks_pos_old;
        update_cluster_disks_pos(:,2:3) = update_coord;
        
        % check overlap for disks of clusters with other clsters and disks
        [condition_cluster] = overlap_check_clusters(ids, update_cluster_disks_pos, disks_stat, dia_disk,neighborList, XX, YY);
        
        if mean(condition_cluster) == 0
           [r,c] = size(update_cluster_disks_pos);
           for kk = 1:r
               iii = update_cluster_disks_pos(kk,1);
               disks_stat(iii,2:3) = update_cluster_disks_pos(kk,2:3);
           end
        end
        
end



