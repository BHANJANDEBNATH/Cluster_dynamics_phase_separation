clear
clc

% Post processing
% Cluster size distribution

size_dstr_samples = [];
pr_dstr_sample = [];

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

    [Ar_a, Pr_a] = size_dstr(disks_stat_time, clusters_stat_time, No_disks);
    size_dstr_samples = [size_dstr_samples Ar_a];
    pr_dstr_sample = [pr_dstr_sample Pr_a];
end

% Mean and errorbar calculation over independent runs

size_mean = zeros(length(size_dstr_samples),1);
size_er = zeros(length(size_dstr_samples),1);

prob_mean = zeros(length(pr_dstr_sample),1);
prob_er = zeros(length(pr_dstr_sample),1);

for ii = 1:length(size_dstr_samples)
    size_mean(ii,1) = mean(size_dstr_samples(ii,1:4));
    size_er(ii,1) = std(size_dstr_samples(ii,1:4));
     
    prob_mean(ii,1) = mean(pr_dstr_sample(ii,1:4));
    prob_er(ii,1) = std(pr_dstr_sample(ii,1:4));
end

zval1 = 1;
zval2 = 2.776;
size_er = zval1 * size_er/4^(0.5);
prob_er = zval2 * prob_er/4^(0.5);


x = size_mean;
xer = size_er;

y = prob_mean;
yer = prob_er;

errorbarxy(x,y,xer,yer,'Color','r','LineStyle','none','Marker','s','LineWidth',2,'MarkerSize',12);

hold on


% % Shaded error area
% x_fill = [x; flipud(x)];
% y_fill = [y+err; flipud(y-err)];
% fill(x_fill, y_fill, [0.5 0.5 1], ...
%      'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
% hold on;
% 
% plot(x, y, 'bo', 'LineWidth', 2)


function [Ar_a, Pr_a] = size_dstr(disks_stat_time, clusters_stat_time, No_disks)

t = length(disks_stat_time);
timeav_samples = 200;

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

ARR = cell(timeav_samples, t/timeav_samples);
ARR_scaled = cell(timeav_samples, t/timeav_samples);

for i = 1:1:(t/timeav_samples)
    for j = 1:1:timeav_samples
        m = (i - 1) * timeav_samples + j;
        A = disks_stat_time{m,2};
        clusters = clusters_stat_time{m,2};
        s = (1:1:No_disks)';
        Ns = zeros(length(s),1);
        for kk = 1:length(clusters)
            cluster_ids = clusters{kk};
            l = length(cluster_ids);
            Ns(l,1)  = Ns(l,1) + 1;
        end
        no_clusters = length(clusters);
        Ns_new = [];
        s_new = [];
        for jj = 1:length(Ns)
            if Ns(jj,1) ~= 0
               Ns_new = [Ns_new; Ns(jj,1)];
               s_new = [s_new; s(jj,1)];
            end
        end
        Ns_dstr = Ns_new/no_clusters;
        Ns_dstr_cumu = Ns_dstr;
        for kkk = 2:length(Ns_dstr)
            Ns_dstr_cumu(kkk,1) = Ns_dstr_cumu(kkk-1,1) + Ns_dstr(kkk,1);
        end

        stat = [s_new Ns_dstr_cumu];
        stat_scaled = [s_new/Vicsek_cluster_size(m,1) Ns_dstr_cumu]; 

        ARR{j,i} = stat;
        ARR_scaled{j,i} = stat_scaled;
    end
end

ARR_arrange = cell(t/timeav_samples, 1);
ARR_scaled_arrange = cell(t/timeav_samples, 1);
for i = 1:1:(t/timeav_samples)
    st = [];
    st_scaled = [];
    for j = 1:timeav_samples
        a = ARR{j,i};
        b = ARR_scaled{j,i};
        st = [st; a];
        st_scaled = [st_scaled; b];
    end
    ARR_arrange{i,1} = st;
    ARR_scaled_arrange{i,1} = st_scaled;
end

TimeAv = cell(t/timeav_samples, 1);

timeav_pr_y = [];
timeav_pr_y_er = [];
timeav_ar_x = [];
timeav_ar_x_er = [];

for i = 1:1:(t/timeav_samples)
    xy = ARR_arrange{i,1};
    x = xy(:,1);        % total area
    y = 1 - xy(:,2);    % (1 - Cumulative Pr)
     
    xy_scaled = ARR_scaled_arrange{i,1};
    x_scaled = xy_scaled(:,1);      % total area
    y_scaled = 1 - xy_scaled(:,2);  % (1 - Cumulative Pr)
        
    binwidth = 0.1;
    pr_min = 0;
    pr_max = 1;
    n_bins = (pr_max - pr_min)/binwidth;
      
    xav = zeros(n_bins,1);
    xav_er = zeros(n_bins,1);
    yav = zeros(n_bins,1);
    yav_er = zeros(n_bins,1);
       
    xav_scaled = zeros(n_bins,1);
    xav_scaled_er = zeros(n_bins,1);
    yav_scaled = zeros(n_bins,1);
    yav_scaled_er = zeros(n_bins,1);
      
    for jj = 1:n_bins
        count = 0;
        xx = [];
        xx_scaled = [];
        yy = [];
        yy_scaled = [];
        for kk = 1:length(y)
            if (y(kk,1) >= (jj - 1) * binwidth) && (y(kk,1) < jj * binwidth)
                count = count + 1;
                xx = [xx; x(kk,1)];
                xx_scaled = [xx_scaled; x_scaled(kk,1)];
                yy = [yy; y(kk,1)];
                yy_scaled = [yy_scaled; y_scaled(kk,1)];
            end
        end
        xav(jj,1) = mean(xx);
        xav_er(jj,1) = std(xx);
        yav(jj,1) = mean(yy);
        yav_er(jj,1) = std(yy);

        xav_scaled(jj,1) = mean(xx_scaled);
        xav_scaled_er(jj,1) = std(xx_scaled);
        yav_scaled(jj,1) = mean(yy_scaled);
        yav_scaled_er(jj,1) = std(yy_scaled);
    end

    f = 1; 
    Stat = [xav (xav_er*f) yav (yav_er*f) xav_scaled (xav_scaled_er*f) yav_scaled (yav_scaled_er*f)];
    TimeAv{i,1} = Stat;
end

time_seq = 10;   % select the time point to extract distribution at that time point
a = TimeAv{time_seq,1};

Ar_a = a(:,5);   % size
Pr_a = a(:,7);   % probability value

end
