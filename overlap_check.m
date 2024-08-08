

% overlap check for disks
function [condition] = overlap_check(id, new_pos_s, dia_disk, sc)

condition = [];

xs = new_pos_s(1);
ys = new_pos_s(2);

% wih other disks
for j  = 1:length(sc)
    if id ~= j
        dis = ((xs - sc(j,1))^2 + (ys - sc(j,2))^2)^0.5;
            if dis > dia_disk
                condition = [condition; 0];
            else
                condition = [condition; 1];
            end
     end
end


end
