

% overlap check for disks
function [condition] = overlap_check_isolated_disks(id, new_pos_s, dia_disk, sc)

condition = [];

xs = new_pos_s(1);
ys = new_pos_s(2);

[r c] = size(sc);

% overlap wih other disks
for j  = 1:r
    if id ~= sc(j,1)
        x2 = sc(j,2);
        y2 = sc(j,3);
        
        % distance
        dis = ((xs - x2)^2 + (ys - y2)^2)^0.5;
            if dis > dia_disk
                condition = [condition; 0];
            else
                condition = [condition; 1];
            end
     end
end


end
