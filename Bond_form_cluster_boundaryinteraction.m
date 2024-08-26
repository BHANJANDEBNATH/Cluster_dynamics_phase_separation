
function[disks_stat] = Bond_form_cluster_boundaryinteraction(No_disks,diskid_out,disks_stat_out,disks_stat,dia_disk, dis_factor, bondpos,bond_form_prob, steps_lifetimes, curr_time, ncols, ncol_bonds,ncol_bond_lifetime,ncol_bond_formtime)


if isempty(diskid_out) ~= 1
    for ii = 1:length(diskid_out)
        id1 = diskid_out(ii,1);
        xo1 = disks_stat_out(id1,2);
        yo1 = disks_stat_out(id1,3);
        for jj = 1:No_disks
            id2 = disks_stat_out(jj,1);
            if id1 ~= id2
               xo2 = disks_stat_out(id2,2);
               yo2 = disks_stat_out(id2,3);
               dis_out = ((xo1 - xo2)^2 + (yo1 - yo2)^2)^0.5;

               % find if any out of 4 postions is empty for bond formation
               mod1 = find(disks_stat_out(id1, bondpos(1,1):bondpos(1,length(bondpos))) == 0);
               mod2 = find(disks_stat_out(id2, bondpos(1,1):bondpos(1,length(bondpos))) == 0);
               % find if any disk id appears repeatedly in bonds
               pod1 = find(disks_stat_out(id1, bondpos(1,1):bondpos(1,length(bondpos))) == id2);
               pod2 = find(disks_stat_out(id2, bondpos(1,1):bondpos(1,length(bondpos))) == id1);
               pr = rand(1,1);   % generate 1 random number between 0 and 1
               if (dis_out > dia_disk && dis_out < dia_disk + dia_disk/dis_factor) && pr < bond_form_prob && isempty(mod1) ~= 1 && isempty(pod1) == 1 && isempty(mod2) ~= 1 && isempty(pod2) == 1
                  % randomly select any one empty position and fill with id of new bonded disk 
                  nod1 = numel(mod1);
                  pod1 = randi(nod1);
                  kkod1 = (ncols - (ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime)) + mod1(pod1);
                  kkkod1 = (ncols - (ncol_bond_lifetime + ncol_bond_formtime)) + mod1(pod1);
                  kkkkod1 = (ncols - (ncol_bond_formtime)) + mod1(pod1);
                  tt = randi(length(steps_lifetimes));
                  disks_stat(id1,kkod1) = id2;
                  disks_stat(id1,kkkod1) = steps_lifetimes(tt,1);
                  disks_stat(id1,kkkkod1) = curr_time;

                  nod2 = numel(mod2);
                  pod2 = randi(nod2);
                  kkod2 = (ncols - (ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime)) + mod2(pod2);
                  kkkod2 = (ncols - (ncol_bond_lifetime + ncol_bond_formtime)) + mod2(pod2);
                  kkkkod2 = (ncols - (ncol_bond_formtime)) + mod2(pod2);
                  disks_stat(id2,kkod2) = id1;
                  disks_stat(id2,kkkod2) = steps_lifetimes(tt,1);
                  disks_stat(id2,kkkkod2) = curr_time;

                  
               end
            end
        end
    end
end


end







