

function [update_cluster_disks_pos] = translate_rotate_cluster(cluster_disks_pos_old, ds, theta_rad)

xxx = cluster_disks_pos_old(:,2)';
yyy = cluster_disks_pos_old(:,3)';

cx = mean(xxx);
cy = mean(yyy);

xt = xxx + ds(1,1);
yt = yyy + ds(1,2);

xnew = [];
ynew = [];

rotate = @(point, angle_rad, center) center + ((point - center) * [cos(angle_rad), -sin(angle_rad); sin(angle_rad), cos(angle_rad)]);

for i = 1:length(xt)
    a = rotate([xt(1,i), yt(1,i)], theta_rad, [cx, cy]);
    xnew = [xnew; a(1)];
    ynew = [ynew; a(2)];
end

update_cluster_disks_pos = cluster_disks_pos_old;
update_cluster_disks_pos(:,2) = xnew;
update_cluster_disks_pos(:,3) = ynew;

end