

% initialization after brownian movement for a few time steps with out any
% interactions

function [updated_coordinates] = initialization_BD_nointeraction(disk_coordinates, XX, YY, dia_disk, tau, D_free_disk, steps, No_disks)

sc = disk_coordinates;

for i = 1:steps

    k = sqrt(2 * 1 * D_free_disk * tau); 
    % updating position
    for j = 1:No_disks
        id = j;
        old_pos_s = sc(id,:);
        ds = k * randn(1,2);
        new_pos_s = old_pos_s + ds;
        
        % periodic BC in x and y directions 
        [new_pos_s] = BoundaryCondition_disks(new_pos_s, XX, YY);

        % checking new_pos_s overlapping with other disks
        [condition] = overlap_check(id, new_pos_s, dia_disk, sc);

        if mean(condition) == 0   % no overlap and update position
            sc(id,:) = new_pos_s;
        end
    end

end

updated_coordinates = sc;



end