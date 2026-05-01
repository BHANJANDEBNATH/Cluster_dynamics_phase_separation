clear
clc


% Post processing
% Neighbor exchanges

Neighbor_exchange_sample = [];
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
    filepath = sprintf('./set%d/data_sim_nb4_af005_N25C1_30000000.mat', i);
    
    load(filepath);
    [c_ex_av, time_av] = neighborexchanges(disks_stat_time,XX);
    Neighbor_exchange_sample = [Neighbor_exchange_sample c_ex_av];
end

% Mean and errorbar calculation over independent runs

NE_mean = zeros(length(Neighbor_exchange_sample),1);
NE_er = zeros(length(Neighbor_exchange_sample),1);

for ii = 1:length(Neighbor_exchange_sample)
    NE_mean(ii,1) = mean(Neighbor_exchange_sample(ii,1:4));
    NE_er(ii,1) = std(Neighbor_exchange_sample(ii,1:4));
end

NE_er = 2.776 * NE_er/4^(0.5);
%V_size_er = V_size_er/4^(0.5);  % standard error

x = time_av;
y = NE_mean;
err = NE_er;

%errorbar(x, y, err,'LineWidth',2)
%hold on

% % % Shaded error area
x_fill = [x; flipud(x)];
y_fill = [y+err; flipud(y-err)];
fill(x_fill, y_fill, [0.5 0.5 1], ...
     'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;

plot(x, y, 'bo', 'LineWidth', 2)







function [c_ex_av, time_av] = neighborexchanges(disks_stat_time,XX)

t = length(disks_stat_time);
L = XX;   
tau = 10^(-5);
time = zeros(t,1);

Neighbor_exchange_av = [];

numbers_exchanged_time = zeros(t,1);

for i = 1:t-1

    A = disks_stat_time{i,2};
    B = disks_stat_time{i+1,2};

    time(i,1) = disks_stat_time{i,1};

    count = 0;
    for j = 1:length(A)
        aa = A(j,4:6);
        bb = B(j,4:6);

        if mean(aa) ~= 0
            aa(aa == 0) = [];
            bb(bb == 0) = [];
            l_aa = length(aa);
            l_bb = length(bb); 
            if l_aa > l_bb
                isSubset = all(ismember(bb, aa));
            else
                isSubset = all(ismember(aa, bb));
            end
            numericValue = double(isSubset);
            if numericValue == 0
                count = count + 1;
            end
        end
    end
    numbers_exchanged_time(i,1) = count; 
end

cumulative_numbers_exchanged = zeros(t,1);

for i = 2:t-1
    a = numbers_exchanged_time(i-1,1);
    cumulative_numbers_exchanged(i,1) = cumulative_numbers_exchanged(i-1,1) + a; 
end

time  = time * 10^(-5);

timeav_samples = 200;
c_ex_av = zeros(t/timeav_samples,1);
c_ex_av_er = zeros(t/timeav_samples,1);
time_av = zeros(t/timeav_samples,1);

for i = 1:1:(t/timeav_samples)
    vv = [];
    tt = [];
    for j = 1:timeav_samples
        m = (i - 1) * timeav_samples + j;
        vv = [vv; cumulative_numbers_exchanged(m,1)];
        tt = [tt; time(m,1)];
    end
    c_ex_av(i,1) = mean(vv);
    c_ex_av_er(i,1) = std(vv);
    time_av(i,1) = mean(tt);
end


end

