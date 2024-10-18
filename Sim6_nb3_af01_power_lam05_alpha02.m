clear
clc
tic;

% Main file 2D simulation of disks

%% simulation box dimensions and other parameters
% XX = 275;           % Height in nm 
% YY = 275;           % Width  in nm 

No_disks = 100;
dia_disk = 10;        % in nm
% area_fraction = ((pi/4) * dia_disk^2 * No_disks) / (XX * YY);

area_fraction = 0.1;      % change area fraction

XX = ((No_disks * pi/4 * dia_disk^2)/area_fraction)^0.5;
XX = round(XX);
YY = XX;
number_bonds = 3;   % number of bonds

%% initial configuration generation
[disk_id, disk_coordinates] = Generate_disk(XX, YY, No_disks, dia_disk);
% Initialization: run BD for a few time steps with out interaction and use the updated
steps = 1000;
tau = 10^(-5);        % in seconds
D_free_disk = 10;   % diffusivity in nm^2/s 
[updated_coordinates] = initialization_BD_nointeraction(disk_coordinates, XX, YY, dia_disk, tau, D_free_disk, steps, No_disks);
% neighborlist initialization
particles = disk_coordinates;
boxSize = [XX, YY];
cutoff = 3 * dia_disk;
[neighborList] = NeighborListCreation(particles, boxSize, cutoff);

%% BD and bond dynamics for disks and clusters
% Bond details ---> allow maximum 4 bonds per disk
% col 1 --> disk id
% col 2, col 3 --> disk x, y  at time t
% col 4, col 5, col 6, col 7 --> ids of disks (maximum 4) that can be bonded with disk id*
% zero represents no bond 

ncol_id = 1;
ncol_pos = 2;      % for two dimensions
ncol_bonds = number_bonds;    % maximum number of bonds
ncol_bond_lifetime = ncol_bonds;   % columns to save bond lifetimes of each bond
ncol_bond_formtime = ncol_bonds;   % columns to save bond formtimes of each bond

ncols = ncol_id + ncol_pos + ncol_bonds + ncol_bond_lifetime + ncol_bond_formtime;

bondpos1 = ncol_id + ncol_pos + 1;
bondpos2 = bondpos1  + ncol_bonds; 
bondpos = bondpos1:1:bondpos2-1;  % column ids to store bonds for any disk

bond_lifetime1 = ncol_id + ncol_pos + ncol_bonds + 1;
bond_lifetime2 = bond_lifetime1  + ncol_bond_lifetime;
bond_lifetime = bond_lifetime1:1:bond_lifetime2-1;   % column ids to store bond lifetime

bond_formtime1 = ncol_id + ncol_pos + ncol_bonds + ncol_bond_lifetime + 1;
bond_formtime2 = bond_formtime1  + ncol_bond_formtime;
bond_formtime =  bond_formtime1:1:bond_formtime2-1;   % column ids to store bond form time

disks_stat = zeros(No_disks, ncols);
disks_stat(:,1) = disk_id;
disks_stat(:,2) = updated_coordinates(:,1);
disks_stat(:,3) = updated_coordinates(:,2);

%% Bond features and lifetimes
lambda = 0.5;          % 0.1 to 2
alpha = 0.2;           % range 0 to 1.5
num_samples = 20000;   % Number of lifetimes to generate
[lifetimes_power, lifetimes_exponential] = distribution_bond_lifetime(alpha, lambda, num_samples);

simulation_timestep = 10^(-5);            % in seconds 

% power-law
steps_lifetimes = lifetimes_power/simulation_timestep;
steps_lifetimes = round(steps_lifetimes);

% exponential
% steps_lifetimes = lifetimes_exponential/simulation_timestep;
% steps_lifetimes = round(steps_lifetimes);

%% Main: Run dynamics
tic_loop = 10^3;

tau = simulation_timestep;    % simulation time step in seconds
D_free_disk = 10^4;           % diffusivity in nm^2/s 
No_timesteps = 60 * 10^5;      
del_t_sampling = 100;

disks_stat_time = cell(No_timesteps/del_t_sampling,2);
clusters_stat_time = cell(No_timesteps/del_t_sampling,2);

for i = 1:No_timesteps
  
    curr_time = i;
    bond_form_prob = 0.5;  % bond formation probability
    dis_factor = 20;       % distance criteria for bond formation
    
    % generate neighborlist after every n time steps to reduce computation
    n_list_steps = 20;
    cutoff_factor = 3; 
    if mod(i,n_list_steps) == 0
        particles = disks_stat(:,2:3);
        boxSize = [XX YY];
        cutoff = cutoff_factor * dia_disk;
        [neighborList] = NeighborListCreation(particles, boxSize, cutoff);
    end
    
    % BD simulation of isolated disks, clusters, and bond formation-breakage dynamics
    [updated_disks_stat, clusters_ids] = BD_disks_bond_dynamics(No_disks, disks_stat, XX, YY, dia_disk, tau, D_free_disk, bond_form_prob, bondpos, bond_lifetime, bond_formtime, steps_lifetimes, curr_time, ncols,ncol_bonds,ncol_bond_lifetime,ncol_bond_formtime,dis_factor,neighborList);
    disks_stat = updated_disks_stat;
    
    % data dumping
    if mod(i,del_t_sampling) == 0
        disks_stat_time{i/del_t_sampling,1} = i;
        disks_stat_time{i/del_t_sampling,2} = disks_stat;
   
        clusters_stat_time{i/del_t_sampling,1} = i;
        clusters_stat_time{i/del_t_sampling,2} = clusters_ids;
    end
    disks_coord_plot = disks_stat;
    
    % checking time elapsed for loop 
    if mod(i,tic_loop) == 0
        disp(['No. time steps elapsed:' num2str(i/tic_loop)]);
    end
    
    % save workspace after every 10 seconds interval
    if mod(i, (10^4 * del_t_sampling)) == 0
        outputname = ['Sim6_nb3_af01_power_lam05_alpha02_' num2str(i) '.mat'];     % change the name of file
        save(outputname)
    end
     
end


%% plotting (orange isolated disks) (purple bonded disks)
% figure(1)
% Plot_disks_colors(No_disks, dia_disk, disks_coord_plot, XX, YY, bondpos)
% 
% % plotting with periodic boxes
% figure(2)
% Plot_disks_colors_with_images(No_disks, dia_disk, disks_coord_plot, XX, YY, bondpos)

elapsedTime = toc;
fprintf('Elapsed time: %.2f seconds\n', elapsedTime);
