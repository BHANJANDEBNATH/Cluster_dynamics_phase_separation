

% Brownian dynamics 
% position update for isolated disks

function[disks_stat] = BD_isolated_disks(ii,disks_stat,D_free_disk,tau,XX,YY,dia_disk)

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