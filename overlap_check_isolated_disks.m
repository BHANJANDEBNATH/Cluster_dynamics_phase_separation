

% overlap check for disks
function [condition] = overlap_check_isolated_disks(id, new_pos_s, dia_disk, sc, XX, YY)

condition = [];

xs = new_pos_s(1);
ys = new_pos_s(2);

% overlap wih other disks
for j  = 1:length(sc)
    if id ~= j
        x2 = sc(j,1);
        y2 = sc(j,2);
        
        % checking whether other disks are inside simulation box or not
        % in x direction
        if x2 <= 0
           x2 = x2 + XX; 
        end
        if x2 > XX
           x2 = x2 - XX;
        end
        % in Y direction
        if y2 <= 0
           y2 = y2 + YY; 
        end
        if y2 > YY
           y2 = y2 - YY;
        end
        
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
