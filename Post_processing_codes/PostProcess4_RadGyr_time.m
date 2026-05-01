clear
clc


% Post processing
% Vicsek cluster size with time 

Rad_gyr_av_samples = [];
for i = 1:4
    % number binding sites 3; average time 130s
    %filepath = sprintf('./set%d/data_sim_nb3_af005_N1C00038_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb3_af005_N5C003_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb3_af005_N25C0047_30000000.mat', i);

    % number binding sites 3; average time 15s
    %filepath = sprintf('./set%d/data_sim_nb3_af005_N1C0034_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb3_af005_N5C03_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb3_af005_N25C1_30000000.mat', i);

    % number binding sites 4; average time 130s
    %filepath = sprintf('./set%d/data_sim_nb4_af005_N1C00038_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb4_af005_N5C003_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb4_af005_N25C0047_30000000.mat', i);

    % number binding sites 4; average time 15s
    %filepath = sprintf('./set%d/data_sim_nb4_af005_N1C0034_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb4_af005_N5C03_30000000.mat', i);
    %filepath = sprintf('./set%d/data_sim_nb4_af005_N25C1_30000000.mat', i);
    
    load(filepath);

    [Rad_Gyr_av, time_av] = Rad_Gyr(disks_stat_time,XX,clusters_stat_time,dia_disk);

    Rad_Gyr_av = Rad_Gyr_av.^0.5;
    Rad_gyr_av_samples = [Rad_gyr_av_samples Rad_Gyr_av];
end

% Mean and errorbar calculation over independent runs

Rad_Gyr_mean = zeros(length(Rad_gyr_av_samples),1);
Rad_Gyr_er = zeros(length(Rad_gyr_av_samples),1);

for ii = 1:length(Rad_gyr_av_samples)
    Rad_Gyr_mean(ii,1) = mean(Rad_gyr_av_samples(ii,1:4));
    Rad_Gyr_er(ii,1) = std(Rad_gyr_av_samples(ii,1:4));
end

Rad_Gyr_er = 2.776 * Rad_Gyr_er/4^(0.5);
%V_size_er = V_size_er/4^(0.5);  % standard error

x = time_av;
y = Rad_Gyr_mean;
err = Rad_Gyr_er;

% Shaded error area
x_fill = [x; flipud(x)];
y_fill = [y+err; flipud(y-err)];
fill(x_fill, y_fill, [0.5 0.5 1], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;

plot(x, y, 'bo', 'LineWidth', 2)




function [Rad_Gyr_av, time_av] = Rad_Gyr(disks_stat_time,XX,clusters_stat_time,dia_disk)

t = length(disks_stat_time);
YY = XX;

tau = 10^(-5);
time = zeros(t,1);
rad_gyr_avg = zeros(t,1);

Mass_wt_Rad_gyr = [];

for i = 1:t
    time(i,1) = disks_stat_time{i,1};
    disks_stat = disks_stat_time{i,2};
    clusters_ids = clusters_stat_time{i,2};
    mass = [];
    Rad_gyr = [];
    for mm = 1:length(clusters_ids)
        ids = clusters_ids{mm,1};
        mmm = length(ids);
        cluster_disks_pos = [];
        if isempty(ids) == 0 && numel(ids) > 1
           for k = 1:length(ids)
               c = [ids(k,1) disks_stat(ids(k,1),2:3)]; 
               cluster_disks_pos = [cluster_disks_pos; c];  % ids, x, y of disks in cluster
           end
           positions = cluster_disks_pos(:,2:3);
           [rad_gyr] = radius_gyration_pbc(ids,positions,XX,YY,dia_disk);
        end 
        Rad_gyr = [Rad_gyr; rad_gyr];
        mass = [mass; mmm];
    end

    m_wt_Rad_gyr = sum(mass .* Rad_gyr.^2)/sum(mass);
    Mass_wt_Rad_gyr = [Mass_wt_Rad_gyr; m_wt_Rad_gyr];
end
    
time = time * tau;

timeav_samples = 200;
Rad_Gyr_av = zeros(t/timeav_samples,1);
time_av = zeros(t/timeav_samples,1);

for i = 1:1:(t/timeav_samples)
    vv = [];
    tt = [];
    for j = 1:timeav_samples
        m = (i - 1) * timeav_samples + j;
        vv = [vv; Mass_wt_Rad_gyr(m,1)];
        tt = [tt; time(m,1)];
    end
    Rad_Gyr_av(i,1) = mean(vv);
    time_av(i,1) = mean(tt);
end


end
