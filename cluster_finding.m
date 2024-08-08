

function [cluster] = cluster_finding(ii, disks_stat, bondpos)

id = ii;
cluster = [];
cluster = [cluster; id];

while 1
    
    l1 = length(cluster);
    for i = 1:l1
        mm = disks_stat(cluster(i,1), bondpos(1,1):bondpos(1,length(bondpos)));
        for j = 1:length(mm)
            if mm(1,j) ~= 0
                % append all ids 
                cluster = [cluster; mm(1,j)];
            end
        end
    end
    
    % remove repeated ids 
    cluster = unique(cluster);
    l2 = length(cluster);

    if l1 - l2 == 0
        break;
    end

end
