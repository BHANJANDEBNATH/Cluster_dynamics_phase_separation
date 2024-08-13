
% BD simulation, bond formation-breakage
% update statistics of disks and clusters at each time step

function [updated_disks_stat, clusters_ids] = BD_disks_bond_dynamics(No_disks,disks_stat, XX, YY, dia_disk, tau, D_free_disk, bond_form_prob, bondpos, bond_lifetime, bond_formtime, steps_lifetimes, curr_time, ncols)


%% Brownian Dynamics of isolated disks
for ii = 1:No_disks
    % if a disk is not bonded with any other disks, then update its position using BD
    if mean(disks_stat(ii, bondpos(1,1):bondpos(1,length(bondpos)))) == 0  
        id = ii;
        old_pos_s = disks_stat(id,2:3);

        % Brownian dynamics and position update
        k = sqrt(2 * 1 * D_free_disk * tau); 
        ds = k * randn(1,2);
        new_pos_s = old_pos_s + ds;

        % periodic BC 
        [new_pos_s] = BoundaryCondition_disks(new_pos_s, XX, YY);

        % check for overlap with other disks who are either bonded or non-bonded
        sc =  disks_stat(:,2:3);
        [condition] = overlap_check_isolated_disks(id, new_pos_s, dia_disk, sc, XX, YY);

        if mean(condition) == 0  
            disks_stat(id,2:3) = new_pos_s;
        end
    end
end

%% bond formation among isolated disks and disks of clusters whose coordinates lie inside simulation box
for ii = 1:No_disks
    % bond formation and store ids of bonded disks
    for jj = 1:No_disks
        if ii ~= jj  
            posx_d1 = disks_stat(ii,2);
            posy_d1 = disks_stat(ii,3);
            posx_d2 = disks_stat(jj,2);
            posy_d2 = disks_stat(jj,3);

            % distance between two disks 
            dis = ((posx_d1 - posx_d2)^2 + (posy_d1 - posy_d2)^2)^0.5;
            % find if any out of 4 postions is empty for bond formation
            mm = find(disks_stat(ii, bondpos(1,1):bondpos(1,length(bondpos))) == 0);
            % find if any disk id appears repeatedly in bonds
            pp = find(disks_stat(ii, bondpos(1,1):bondpos(1,length(bondpos))) == jj);
            pr = rand(1,1);   % generate 1 random number between 0 and 1
            dis_factor = 20;
            if (dis > dia_disk && dis < dia_disk + dia_disk/dis_factor) && pr < bond_form_prob && isempty(mm) ~= 1 && isempty(pp) == 1 
                % randomly select any one empty position and fill with id of new bonded disk 
                nn = numel(mm);
                pp = randi(nn);
                kk = 3 + mm(pp);
                kkk = 7 + mm(pp);
                kkkk = 11 + mm(pp);
                disks_stat(ii,kk) = jj;
                tt = randi(length(steps_lifetimes));
                disks_stat(ii,kkk) = steps_lifetimes(tt,1);
                disks_stat(ii,kkkk) = curr_time;
            end
        end
    end
end

%% bond break based on lifetime
for aa = 1:No_disks
    for bb = bondpos(1,1):1:bondpos(1,length(bondpos))
        if disks_stat(aa,bb) ~= 0 
            colid_lt = bb + length(bond_lifetime);
            colid_ft = colid_lt + length(bond_formtime);
            lt = disks_stat(aa,colid_lt);   % lifetime of bonded disk
            ft = disks_stat(aa,colid_ft);   % time at which bond formed 
            tt = (curr_time - ft);          % time elapsed after bond formation
            if tt == lt
                disks_stat(aa,bb) = 0;
                disks_stat(aa,colid_lt)  = 0;
                disks_stat(aa,colid_ft)  = 0;
            end
        end
    end
end

%% finding clusters
clusters_ids = cell(No_disks,1);
for ii = 1:No_disks
    % if a disk is bonded, then find the associated cluster
    if mean(disks_stat(ii, bondpos(1,1):bondpos(1,length(bondpos)))) ~= 0 
        % cluster finding operation
        [cluster] = cluster_finding(ii, disks_stat, bondpos);
        clusters_ids{ii} = cluster;
    end
end
% remove empty and repeating cells of clusters
cc = clusters_ids;
ccnew = cc(~cellfun(@isempty, cc));
F = ccnew;
for ii = numel(F):-1:1
    for jj = 1:ii-1
        if isequal(F{ii},F{jj})
           F(ii) = [];
           break
        end
     end
end
% final sets of clusters
clusters_ids = F;

%% move cluster based on its center of mass and 
% update positions of disks associated with the clusters
for m = 1:length(clusters_ids)    
    ids = clusters_ids{m,1};
    cluster_disks_pos_old = [];
    update_cluster_disks_pos = [];

    if isempty(ids) == 0 && numel(ids) > 1

        for k = 1:length(ids)
            c_old = [ids(k,1) disks_stat(ids(k,1),2:3)]; 
            cluster_disks_pos_old = [cluster_disks_pos_old; c_old];  % ids, x, y of disks in cluster
        end

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
          
end

%% update final positions of all disks
updated_disks_stat = disks_stat;

end




