
% Bond formation among isolated disks and disks belonged to clusters 
% whose coordinates lie inside simulation box
% no bond formation interaction across boundaries
% ISB --- inside simulation box

function[disks_stat] = Bond_form_ISB(id_d1,posx_d1,posy_d1,id_d2,posx_d2,posy_d2, disks_stat,bondpos,ncols,ncol_bonds,ncol_bond_lifetime,ncol_bond_formtime,steps_lifetimes,curr_time, dia_disk,dis_factor,bond_form_prob)


            % distance between two disks 
            dis = ((posx_d1 - posx_d2)^2 + (posy_d1 - posy_d2)^2)^0.5;
            % find if any out of 4 postions is empty for bond formation
            mm_d1 = find(disks_stat(id_d1, bondpos(1,1):bondpos(1,length(bondpos))) == 0);
            mm_d2 = find(disks_stat(id_d2, bondpos(1,1):bondpos(1,length(bondpos))) == 0);
            % find if any disk id appears repeatedly in bonds
            pp_d1 = find(disks_stat(id_d1, bondpos(1,1):bondpos(1,length(bondpos))) == id_d2);
            pp_d2 = find(disks_stat(id_d2, bondpos(1,1):bondpos(1,length(bondpos))) == id_d1);
            pr = rand(1,1);   % generate 1 random number between 0 and 1
            if (dis > dia_disk && dis < dia_disk + dia_disk/dis_factor) && pr < bond_form_prob && isempty(mm_d1) ~= 1 && isempty(pp_d1) == 1 && isempty(mm_d2) ~= 1 && isempty(pp_d2) == 1 
                % randomly select any one empty position and fill with id of new bonded disk 
                nn_d1 = numel(mm_d1);
                pp_d1 = randi(nn_d1);
                kk_d1 = (ncols - (ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime)) + mm_d1(pp_d1);
                kkk_d1 = (ncols - (ncol_bond_lifetime + ncol_bond_formtime)) + mm_d1(pp_d1);
                kkkk_d1 = (ncols - (ncol_bond_formtime)) + mm_d1(pp_d1);
                disks_stat(id_d1,kk_d1) = id_d2;
                tt = randi(length(steps_lifetimes));
                disks_stat(id_d1,kkk_d1) = steps_lifetimes(tt,1);
                disks_stat(id_d1,kkkk_d1) = curr_time;

                nn_d2 = numel(mm_d2);
                pp_d2 = randi(nn_d2);
                kk_d2 = (ncols - (ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime)) + mm_d2(pp_d2);
                kkk_d2 = (ncols - (ncol_bond_lifetime + ncol_bond_formtime)) + mm_d2(pp_d2);
                kkkk_d2 = (ncols - (ncol_bond_formtime)) + mm_d2(pp_d2);
                disks_stat(id_d2,kk_d2) = id_d1;
                disks_stat(id_d2,kkk_d2) = steps_lifetimes(tt,1);
                disks_stat(id_d2,kkkk_d2) = curr_time;
            end



end






