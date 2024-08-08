clear
clc

% Main 2D simulation of disks

%% dimensions 
XX = 100;    % Height in nm 
YY = 100;    % Width  in nm 
No_disks = 150;
dia_disk = 5;   % in nm
area_fraction = ((pi/4) * dia_disk^2 * No_disks) / (XX * YY); 

%% initial configuration generation
[disk_id,disk_coordinates] = Generate_disk(XX,YY, No_disks, dia_disk);
%figure
%Plot_disks(No_disks, dia_disk, disk_coordinates)

%% Initialization: run BD for a few time steps with out interaction and use the updated
% coordinates as new initial configuration
tau = 10^(-3);        % in s
D_free_disk = 10^2;   % diffusivity in nm^2/s (diffusivity 10^(-14) m^2/s)

steps = 10^3;
[updated_coordinates] = initialization_BD_nointeraction(disk_coordinates, XX, YY, dia_disk, tau, D_free_disk, steps, No_disks);
%figure
%Plot_disks(No_disks, dia_disk, updated_coordinates)

%% BD and bond dynamics for disks and clusters
% Bond details ---> allow maximum 4 bonds per disk
% col 1 --> disk id*
% col 2, col 3 --> disk x, y  at time t
% col 4, col 5, col 6, col 7 --> ids of disks (maximum 4) that can be bonded with disk id*
% zero represents no bond 

ncol_id = 1;
ncol_pos = 2;      % for two dimensions
ncol_bonds = 4;    % allowing maximum 4 bonds
ncol_bond_lifetime = ncol_bonds;
ncol_bond_formtime = ncol_bonds;

ncols = ncol_id + ncol_pos + ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime;
bondpos = linspace(ncols - ncol_id - ncol_pos - ncol_bond_lifetime - ncol_bond_formtime, ncols - ncol_bond_lifetime - ncol_bond_formtime, ncol_bonds);  % column ids to store bonds for any disk
bond_lifetime = linspace(ncols - ncol_id - ncol_pos - ncol_bond_formtime, ncols - ncol_bond_formtime, ncol_bond_lifetime);  % column ids to store bond lifetime
bond_formtime = linspace(ncols - ncol_id - ncol_pos, ncols, ncol_bond_formtime);  % column id to store bond lifetime

disks_stat = zeros(No_disks, ncols);
disks_stat(:,1) = disk_id;
disks_stat(:,2) = updated_coordinates(:,1);
disks_stat(:,3) = updated_coordinates(:,2);

% Bond lifetimes
alpha = 0.5;          % range 0 1.5
lambda = 0.1;         % 0.1 to 2
num_samples = 1000;  % Number of lifetimes to generate
[lifetimes_power, lifetimes_exponential] = distribution_bond_lifetime(alpha, lambda, num_samples);

%steps_lifetimes = lifetimes_power/tau;
%steps_lifetimes = round(steps_lifetimes);
steps_lifetimes = lifetimes_exponential/tau;
steps_lifetimes = round(steps_lifetimes);

%% Run dynamics
tic = 10^4;
No_timesteps = 1*10^4;
del_t_sampling = 10;
disks_stat_time = cell(No_timesteps/del_t_sampling,2);
clusters_stat_time = cell(No_timesteps/del_t_sampling,2);

for i = 1:No_timesteps

    curr_time = i;
    bond_form_prob = 0.8;  % bond formation probability
    [updated_disks_stat, clusters_ids] = BD_disks_bond_dynamics(No_disks, disks_stat, XX, YY, dia_disk, tau, D_free_disk, bond_form_prob, bondpos, bond_lifetime, bond_formtime, steps_lifetimes, curr_time, ncols);
    disks_stat = updated_disks_stat;
    
    % data sampling
    if mod(i,del_t_sampling) == 0
        disks_stat_time{i/del_t_sampling,1} = i;
        disks_stat_time{i/del_t_sampling,2} = disks_stat;
        
        clusters_stat_time{i/del_t_sampling,1} = i;
        clusters_stat_time{i/del_t_sampling,2} = clusters_ids;
    end
      
    % time check
    if mod(i,tic) == 0
        % tracking how many time steps are elapsed
        disp(['No. time steps elapsed:' num2str(i/tic)]);
    end

    Plot_disks_colors(No_disks, dia_disk, disks_stat, XX, YY, bondpos)

end

%% plotting 
% blue isolated disks
% red bonded disks
% Plot_disks_colors(No_disks, dia_disk, disks_stat, XX, YY, bondpos)













