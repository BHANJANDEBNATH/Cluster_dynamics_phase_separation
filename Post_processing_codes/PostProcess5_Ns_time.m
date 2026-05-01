clear
clc


% Post processing
% Vicsek cluster size with time 

Ns_vicsek_av_samples = [];

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

    [Ns_vicsek_av, time_av] = clustersize(disks_stat_time,XX,clusters_stat_time,No_disks);
    Ns_vicsek_av_samples = [Ns_vicsek_av_samples Ns_vicsek_av];
end

% Mean and errorbar calculation over independent runs

%Ns_vicsek_av_samples = Ns_vicsek_av_samples * XX^2;

Ns_mean = zeros(length(Ns_vicsek_av_samples),1);
Ns_er = zeros(length(Ns_vicsek_av_samples),1);

for ii = 1:length(Ns_vicsek_av_samples)
    Ns_mean(ii,1) = mean(Ns_vicsek_av_samples(ii,1:4));
    Ns_er(ii,1) = std(Ns_vicsek_av_samples(ii,1:4));
end

Ns_er = 2.776 * Ns_er/4^(0.5);
%V_size_er = V_size_er/4^(0.5);  % standard error

x = time_av;
y = Ns_mean;
err = Ns_er;

% Shaded error area
x_fill = [x; flipud(x)];
y_fill = [y+err; flipud(y-err)];
fill(x_fill, y_fill, [0.5 0.5 1], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;

plot(x, y, 'bo-', 'LineWidth', 2)







function [Ns_vicsek_av, time_av] = clustersize(disks_stat_time,XX,clusters_stat_time,No_disks)

t = length(disks_stat_time);
L = XX;   
tau = 10^(-5);
time = zeros(t,1);

Ns_vicsek = [];

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

    Ns_vicsek = [Ns_vicsek; sum(Ns)];
end

time = time * tau;

timeav_samples = 200;
time_av = zeros(t/timeav_samples,1);
Ns_vicsek_av = zeros(t/timeav_samples,1);

for i = 1:1:(t/timeav_samples)
    vv = [];
    tt = [];
    for j = 1:timeav_samples
        m = (i - 1) * timeav_samples + j;
        vv = [vv; Ns_vicsek(m,1)];
        tt = [tt; time(m,1)];
    end
    Ns_vicsek_av(i,1) = mean(vv);
    time_av(i,1) = mean(tt);
end


end

