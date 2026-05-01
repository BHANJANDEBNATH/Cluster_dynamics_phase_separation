
% LAMMPS dump file generation

%load('data_sim_nb3_af01_N1C0034_30000000.mat')
%load('data_sim_nb3_af01_N5C03_30000000.mat')


t = length(disks_stat_time);
L = XX;   
tau = 10^(-5);
time = zeros(t,1);


for i = 1:t
     
    time(i,1) = disks_stat_time{i,1};
    A = disks_stat_time{i,2};

    dump = zeros(length(A),5);
    count = 0;
    
    mm = i*100;
    fid = fopen("BR1." + mm,'w');
    fprintf(fid, 'ITEM: TIMESTEP\n');
    fprintf(fid, '%d\n',i);
    fprintf(fid, 'ITEM: NUMBER OF ATOMS\n');
    fprintf(fid, '%d\n',No_disks);
    fprintf(fid, 'ITEM: BOX BOUNDS pp pp pp\n');
    fprintf(fid, '%f %f\n', [0, L]);
    fprintf(fid, '%f %f\n', [0, L]);
    fprintf(fid, '%f %f\n', [-1, 1]);
    fprintf(fid, 'ITEM: ATOMS id type x y z\n');
    for jj = 1:length(A)
        count = count + 1;
        fprintf(fid, '%d %d %f %f %f\n', [count, count, A(jj,2), A(jj,3)], 0);
    end
    fclose(fid);

end
