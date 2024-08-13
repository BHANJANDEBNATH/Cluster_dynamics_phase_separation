
% overlap check for clusters

function [condition_cluster] = overlap_check_clusters(ids, update_cluster_disks_pos, disks_stat, dia_disk, XX, YY)

condition_cluster = [];

for i = 1:length(ids)
    posx = update_cluster_disks_pos(i,2);
    posy = update_cluster_disks_pos(i,3);
    
    % for overlap check among all disks, it is mandatory that all disks to be
    % inside the simulation box
    % it may happen that though the center of mass of cluster is inside simulation box
    % but some disks corrsponding to that cluster are outside the simulation box 
    % because we are treating each cluster as one large entity
    % note the center of mass of cluster is always inside the simulation box (see overlap_check_clusters.m)
    % that's why, for overlap check for clusters, shifting positions is mandatory for those disks of a cluster 
    % whose coordinates are outside the simulation domain  
    
    % first shift coordiantes of disk inside simulation box
    % if the disk stays outside the box 
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
    
    % second, check overlap wih other disks
    % check whether other disks are inside simulation box or not
    for j = 1:length(disks_stat)
        if ids(i,1) ~= disks_stat(j,1)
            xx2 = disks_stat(j,2);
            yy2 = disks_stat(j,3);

            % checking whether other disks are inside simulation box or not
            % in x direction
            if xx2 <= 0
               xx2 = xx2 + XX; 
            end
            if xx2 > XX
               xx2 = xx2 - XX;
            end
            % in Y direction
            if yy2 <= 0
               yy2 = yy2 + YY; 
            end
            if yy2 > YY
               yy2 = yy2 - YY;
            end
            
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