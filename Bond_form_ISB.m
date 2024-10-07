
% Bond formation among isolated disks and disks belonged to clusters 
% whose coordinates lie inside simulation box
% no bond formation interaction across boundaries
% ISB --- inside simulation box

function[disks_stat] = Bond_form_ISB(id_d1,posx_d1,posy_d1,id_d2,posx_d2,posy_d2, disks_stat,bondpos,ncols,ncol_bonds,ncol_bond_lifetime,ncol_bond_formtime,steps_lifetimes,curr_time, dia_disk,dis_factor,bond_form_prob,XX,YY)

            % distance between two disks
            % consider periodic boundary interactions
            % cosider second disk and create its image-coordinates and measure distances 
            dis1 = ((posx_d1 - posx_d2)^2 + (posy_d1 - posy_d2)^2)^0.5;
            dis_image2 = ((posx_d1 - (posx_d2+XX))^2 + (posy_d1 - posy_d2)^2)^0.5;
            dis_image3 = ((posx_d1 - (posx_d2+XX))^2 + (posy_d1 - (posy_d2+YY))^2)^0.5;
            dis_image4 = ((posx_d1 - (posx_d2+XX))^2 + (posy_d1 - (posy_d2-YY))^2)^0.5;
            dis_image5 = ((posx_d1 - (posx_d2))^2 + (posy_d1 - (posy_d2+YY))^2)^0.5;
            dis_image6 = ((posx_d1 - (posx_d2))^2 + (posy_d1 - (posy_d2-YY))^2)^0.5;
            dis_image7 = ((posx_d1 - (posx_d2-XX))^2 + (posy_d1 - posy_d2)^2)^0.5;
            dis_image8 = ((posx_d1 - (posx_d2-XX))^2 + (posy_d1 - (posy_d2+YY))^2)^0.5;
            dis_image9 = ((posx_d1 - (posx_d2-XX))^2 + (posy_d1 - (posy_d2-YY))^2)^0.5;

            Dis = [dis1; dis_image2; dis_image3; dis_image4; dis_image5; dis_image6; dis_image7; dis_image8; dis_image9];
            dis = min(Dis);
            
            % find if any out of 4 postions is empty for bond formation
            vec1 = disks_stat(id_d1, bondpos(1,1):bondpos(1,length(bondpos)));
            % mm_d1 will store the colums ids of vec1 containing zeros
            mm_d1 = find(vec1 == 0);
            
            vec2 = disks_stat(id_d2, bondpos(1,1):bondpos(1,length(bondpos)));
            mm_d2 = find(vec2 == 0);

            % find if any disk id appears repeatedly in bonds
            nn_d1 = find(vec1 == id_d2);
            nn_d2 = find(vec2 == id_d1);

            % if isempty(mm_d1) = 0, i.e. bond formation allowed because total number of bonds is not exceeded 
            % if isempty(mm_d1) = 1, i.e. bond formaton cannot be allowed
            % if isempty(nn_d1) = 0, i.e. id_d2 is already bonded with id1, hence avoid the loop
            % if isempty(nn_d1) = 1, i.e. id_d2 is not bonded and perform the operation

            pr = rand(1,1);   % generate 1 random number between 0 and 1
            if (dis > dia_disk && dis < dia_disk + dia_disk/dis_factor) && pr < bond_form_prob && isempty(mm_d1) ~= 1 && isempty(nn_d1) == 1 && isempty(mm_d2) ~= 1 && isempty(nn_d2) == 1 
                
                % to randomly select any one empty position and fill with id of new bonded disk 
                nn_d1 = numel(mm_d1);
                pp_d1 = randi(nn_d1); 
                c1 = mm_d1(pp_d1); 

                % column ids
                kk_d1 = (ncols - (ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime)) + c1; % column id to store id of the bonded particle
                kkk_d1 = (ncols - (ncol_bond_lifetime + ncol_bond_formtime)) + c1;  % column id to store bond life time of the corresponding bonded particle
                kkkk_d1 = (ncols - (ncol_bond_formtime)) + c1;    % column id to store bond formation time of the corresponding bonded particle
       
                disks_stat(id_d1,kk_d1) = id_d2;
                tt = randi(length(steps_lifetimes));
                disks_stat(id_d1,kkk_d1) = steps_lifetimes(tt,1);
                disks_stat(id_d1,kkkk_d1) = curr_time;
                
                nn_d2 = numel(mm_d2);
                pp_d2 = randi(nn_d2);
                c2 = mm_d2(pp_d2);

                % column ids
                kk_d2 = (ncols - (ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime)) + c2;
                kkk_d2 = (ncols - (ncol_bond_lifetime + ncol_bond_formtime)) + c2;
                kkkk_d2 = (ncols - (ncol_bond_formtime)) + c2;
                
                disks_stat(id_d2,kk_d2) = id_d1;
                disks_stat(id_d2,kkk_d2) = steps_lifetimes(tt,1);
                disks_stat(id_d2,kkkk_d2) = curr_time;
            end   


end






