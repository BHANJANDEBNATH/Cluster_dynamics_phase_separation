function[disks_coord_plot] = disks_coord_simbox_plot(disks_stat,XX,YY)

disks_coord_plot = disks_stat;

for ii = 1:length(disks_stat)
    id = disks_stat(ii,1);
    posx = disks_stat(ii,2);
    posy = disks_stat(ii,3);

    new_pos_s = [posx posy];
    % use periodic boundary condition and shift coordinates first 
    [new_pos_s] = BoundaryCondition_disks(new_pos_s, XX, YY);
    disks_coord_plot(id,2) = new_pos_s(1);
    disks_coord_plot(id,3) = new_pos_s(2);
end

end