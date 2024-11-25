
% bond break based on lifetime

function[disks_stat] = Bond_break_lifetime(ii,bondpos,disks_stat,bond_lifetime,bond_formtime,curr_time)

    aa = ii;
    mm = bondpos(1,1);
    mmm = bondpos(1,length(bondpos));
    for bb = mm:1:mmm
        if disks_stat(aa,bb) ~= 0 
            colid_lt = bb + length(bond_lifetime);
            colid_ft = colid_lt + length(bond_formtime);
            lt = disks_stat(aa,colid_lt);   % lifetime of bonded disk
            ft = disks_stat(aa,colid_ft);   % time at which bond formed 
            tt = (curr_time - ft);          % time elapsed after bond formation
            if tt == lt
                disks_stat(aa,bb) = 0;
                disks_stat(aa,colid_lt)  = 0;
                disks_stat(aa,colid_ft)  = 0;
            end
        end
    end

end