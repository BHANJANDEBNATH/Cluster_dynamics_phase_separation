% function file generating coordinates of spheres

function [disk_id,disk_coordinates] = Generate_disk(XX,YY, No_disks, dia_disk)

% generating spheres at random positions non overlapping with fibers
disk_coordinates = zeros(No_disks,2);
disk_id = zeros(No_disks,1);

id = 0;
for i = 1:No_disks

    while 1
        xs = (XX - 0) * rand(1);
        ys = (YY - 0) * rand(1);
        
        condition = [];

        for j  = 1:length(disk_coordinates)
            dis = ((xs - disk_coordinates(j,1))^2 + (ys - disk_coordinates(j,2))^2)^0.5;
            if dis >  dia_disk
                condition = [condition; 0];
            else
                condition = [condition; 1];
            end
        end


        if mean(condition) == 0   % if all zero,  then no overlap with other disks
            break;
        end

    end

    disk_coordinates(i,:) = [xs ys];
    id = id + 1;
    disk_id(i,1) = id;

end

end