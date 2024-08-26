

function[diskid_out,disks_stat_out] = disks_clusters_outside(No_disks,disks_stat_out,XX,YY,diskid_out)

for ii = 1:No_disks
    id = disks_stat_out(ii,1);
    posx = disks_stat_out(ii,2);
    posy = disks_stat_out(ii,3);
    if posx <= 0 || posx > 0 || posy <= 0 ||posy > 0
        diskid_out = [diskid_out; id];
    end 
    new_pos_s = [posx posy];
    % use periodic boundary condition and shift coordinates first 
    [new_pos_s] = BoundaryCondition_disks(new_pos_s, XX, YY);
    disks_stat_out(id,2) = new_pos_s(1);
    disks_stat_out(id,3) = new_pos_s(2);
end


end