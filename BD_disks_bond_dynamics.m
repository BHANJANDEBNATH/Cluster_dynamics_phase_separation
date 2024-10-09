
% BD simulation, bond formation-breakage
% update statistics of disks and clusters at each time step

function [updated_disks_stat, clusters_ids] = BD_disks_bond_dynamics(No_disks, disks_stat, XX, YY, dia_disk, tau, D_free_disk, bond_form_prob, bondpos, bond_lifetime, bond_formtime, steps_lifetimes, curr_time, ncols,ncol_bonds,ncol_bond_lifetime,ncol_bond_formtime,dis_factor,neighborList)

%% Brownian Dynamics of isolated disks
for ii = 1:No_disks
    % if a disk is not bonded with any other disks, then update its position using BD
    if mean(disks_stat(ii, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        [disks_stat] = BD_isolated_disks(ii,disks_stat,D_free_disk,tau,XX,YY,dia_disk,neighborList);
    end
end

%% bond formation among isolated disks and disks of clusters
for ii = 1:No_disks
    % bond formation and store ids of bonded disks 
    % perform bond formation for neghboring particles
    neighbor_ids = neighborList{ii,1};
    if isempty(neighbor_ids) ~= 1
      for jj = 1:length(neighbor_ids)
         if ii ~= jj 
            id_d1 = disks_stat(ii,1);
            posx_d1 = disks_stat(ii,2);
            posy_d1 = disks_stat(ii,3);

            n_id = neighbor_ids(jj); 
            id_d2 = disks_stat(n_id,1);
            posx_d2 = disks_stat(n_id,2);
            posy_d2 = disks_stat(n_id,3);
            [disks_stat] = Bond_form_ISB(id_d1,posx_d1,posy_d1,id_d2,posx_d2,posy_d2, disks_stat,bondpos,ncols,ncol_bonds,ncol_bond_lifetime,ncol_bond_formtime,steps_lifetimes,curr_time, dia_disk,dis_factor,bond_form_prob,XX,YY);
         end
      end
    end
end

% bond break based on lifetime
for ii = 1:No_disks
    [disks_stat] = Bond_break_lifetime(ii,bondpos,disks_stat,bond_lifetime,bond_formtime,curr_time);
end

% finding clusters
clusters_ids = cell(No_disks,1);
for ii = 1:No_disks
    % if a disk is bonded, then find the associated cluster
    if mean(disks_stat(ii, bondpos(1,1):bondpos(1,length(bondpos)))) ~= 0 
        % cluster finding operation
        [cluster] = cluster_finding(ii, disks_stat, bondpos);
        clusters_ids{ii} = cluster;
    end
end
% Arranging --- remove empty and repeating cells of clusters
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

%% move cluster based on its center of mass 
% update positions of disks associated with the clusters
for m = 1:length(clusters_ids)    
    ids = clusters_ids{m,1};
    cluster_disks_pos_old = [];
    if isempty(ids) == 0 && numel(ids) > 1
        for k = 1:length(ids)
            c_old = [ids(k,1) disks_stat(ids(k,1),2:3)]; 
            cluster_disks_pos_old = [cluster_disks_pos_old; c_old];  % ids, x, y of disks in cluster
        end
        [disks_stat] = BD_Clusters_pos_update(ids,cluster_disks_pos_old,disks_stat,D_free_disk,dia_disk,tau,XX,YY,neighborList);        
    end    
end

%% update final positions of all disks
updated_disks_stat = disks_stat;

end




