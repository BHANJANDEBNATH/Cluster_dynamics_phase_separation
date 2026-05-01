clear
clc


% Post processing
% Vicsek cluster size with time 

Vicsek_size_samples = [];
for i = 1:4
    % number binding sites 3; average time 130s
    filepath = sprintf('./set%d/data_sim_nb3_af005_N1C00038_30000000.mat', i);
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

    [Vicsek_cluster_size_av, time_av] = clustersize(disks_stat_time,XX,clusters_stat_time,No_disks);
    Vicsek_size_samples = [Vicsek_size_samples Vicsek_cluster_size_av];
end

% Mean and errorbar calculation over independent runs

V_size_mean = zeros(length(Vicsek_size_samples),1);
V_size_er = zeros(length(Vicsek_size_samples),1);

for ii = 1:length(Vicsek_size_samples)
    V_size_mean(ii,1) = mean(Vicsek_size_samples(ii,1:4));
    V_size_er(ii,1) = std(Vicsek_size_samples(ii,1:4));
end

V_size_er = 2.776 * V_size_er/4^(0.5);
%V_size_er = V_size_er/4^(0.5);  % standard error

x = time_av;
y = V_size_mean;
err = V_size_er;

% Shaded error area
x_fill = [x; flipud(x)];
y_fill = [y+err; flipud(y-err)];
fill(x_fill, y_fill, [0.5 0.5 1], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;

plot(x, y, 'bo', 'LineWidth', 2)







function [Vicsek_cluster_size_av, time_av] = clustersize(disks_stat_time,XX,clusters_stat_time,No_disks)

t = length(disks_stat_time);
L = XX;   
tau = 10^(-5);
time = zeros(t,1);

Vicsek_cluster_size = [];

for i = 1:t
    time(i,1) = disks_stat_time{i,1};
    A = disks_stat_time{i,2};
    clusters = clusters_stat_time{i,2};
    s = (1:1:No_disks)';
    Ns = zeros(length(s),1);

    for k = 1:length(clusters)
        cluster_ids = clusters{k};
        l = length(cluster_ids);
        Ns(l,1)  = Ns(l,1) + 1;
    end

    size = sum(Ns.* s.^2)/sum(Ns.* s);
    Vicsek_cluster_size = [Vicsek_cluster_size; size];
end

time = time * tau;

timeav_samples = 200;
Vicsek_cluster_size_av = zeros(t/timeav_samples,1);
Vicsek_cluster_size_av_er = zeros(t/timeav_samples,1);
time_av = zeros(t/timeav_samples,1);

for i = 1:1:(t/timeav_samples)
    vv = [];
    tt = [];
    for j = 1:timeav_samples
        m = (i - 1) * timeav_samples + j;
        vv = [vv; Vicsek_cluster_size(m,1)];
        tt = [tt; time(m,1)];
    end
    Vicsek_cluster_size_av(i,1) = mean(vv);
    Vicsek_cluster_size_av_er(i,1) = std(vv);
    time_av(i,1) = mean(tt);
end


end

































