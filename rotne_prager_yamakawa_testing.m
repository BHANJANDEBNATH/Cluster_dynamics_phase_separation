function brownian_dynamics_rpy_ewald()
    % Parameters
    N = 100;         % Number of particles
    L = 10.0;        % Size of the simulation box
    alpha = 1.0;     % Splitting parameter for Ewald summation
    dt = 0.01;       % Time step for Brownian Dynamics
    gamma = 1.0;     % Friction coefficient
    num_steps = 1000; % Number of simulation steps
    radius = 0.5;    % Radius of the particles (assumed spherical)

    % Initialize particle positions and velocities
    positions = L * rand(N, 3);
    velocities = randn(N, 3);

    % Main simulation loop
    for step = 1:num_steps
        % Compute pairwise distances
        distances = compute_distances(positions, L);

        % Compute the Ewald summation for long-range hydrodynamic interactions
        ewald_sum = ewald_summation(distances, alpha, L);

        % Compute Rotne-Prager-Yamakawa tensor
        rpy_tensor = rpy_tensor(distances, radius);

        % Update positions and velocities
        [positions, velocities] = update_positions_velocities(positions, velocities, dt, gamma, L, ewald_sum, rpy_tensor);

        % Print progress
        if mod(step, 100) == 0
            fprintf('Step %d: Positions mean = %.2f, Velocities mean = %.2f\n', step, mean(positions(:)), mean(velocities(:)));
        end
    end
end

function distances = compute_distances(positions, L)
    % Compute pairwise distances with periodic boundary conditions
    N = size(positions, 1);
    distances = zeros(N, N);
    for i = 1:N
        for j = 1:N
            if i ~= j
                diff = positions(i, :) - positions(j, :);
                diff = diff - L * round(diff / L); % Apply periodic boundary conditions
                distances(i, j) = sqrt(sum(diff.^2));
            end
        end
    end
end

function ewald_sum = ewald_summation(distances, alpha, L)
    % Compute the Ewald summation for long-range hydrodynamic interactions
    short_range_cutoff = 2.5; % Example cutoff distance for short-range interactions
    short_range = distances;
    short_range(short_range >= short_range_cutoff) = 0;
    short_range_sum = sum(1 ./ short_range(short_range > 0), 'all');
    
    % Reciprocal space summation
    recip_space_sum = 0;
    range = -5:5;
    [kx, ky, kz] = ndgrid(range, range, range);
    kx = kx(:);
    ky = ky(:);
    kz = kz(:);
    for i = 1:length(kx)
        k = [kx(i), ky(i), kz(i)] * 2 * pi / L;
        k2 = sum(k.^2);
        if k2 > 0
            recip_space_sum = recip_space_sum + exp(-k2 / (4 * alpha^2)) / k2;
        end
    end

    ewald_sum = (short_range_sum + recip_space_sum) * (L^3 / (4 * pi));
end

function rpy_tensor = rpy_tensor(distances, radius)
    % Compute the Rotne-Prager-Yamakawa tensor
    N = size(distances, 1);
    rpy_tensor = zeros(N, N, 3, 3);
    for i = 1:N
        for j = 1:N
            if i ~= j
                r = distances(i, j);
                if r > 0
                    % Rotne-Prager-Yamakawa tensor components
                    A = (1 / r) * (1 + 0.5 * radius / r);
                    B = -1 / (6 * r^2);
                    rpy_tensor(i, j, :, :) = A * eye(3) + B * (eye(3) - ones(3, 3));
                end
            end
        end
    end
end

function [positions, velocities] = update_positions_velocities(positions, velocities, dt, gamma, L, ewald_sum, rpy_tensor)
    % Update velocities and positions
    N = size(positions, 1);
    new_velocities = zeros(N, 3);
    
    % Compute hydrodynamic forces and update velocities
    for i = 1:N
        force = zeros(1, 3);
        for j = 1:N
            if i ~= j
                % Rotne-Prager-Yamakawa force contribution
                r = compute_distances(positions([i, j], :), L);
                if r > 0
                    force = force + squeeze(rpy_tensor(i, j, :, :)) * (positions(i, :) - positions(j, :))';
                end
            end
        end
        % Update velocities (Brownian motion)
        new_velocities(i, :) = velocities(i, :) + (force / gamma) * dt;
    end

    % Update positions
    positions = positions + new_velocities * dt;

    % Apply periodic boundary conditions
    positions = mod(positions, L);
    velocities = new_velocities;
end



% for ellipsoidal

% function ellipsoidal_rpy_tensor(distances, a, b, c)
%     % Compute the Rotne-Prager-Yamakawa tensor for ellipsoidal particles
%     % a, b, c: semi-axes of the ellipsoid
% 
%     N = size(distances, 1);
%     rpy_tensor = zeros(N, N, 3, 3);
% 
%     for i = 1:N
%         for j = 1:N
%             if i ~= j
%                 r = distances(i, j);
%                 if r > 0
%                     % Example formula for ellipsoidal particles
%                     % Requires the use of special functions for exact calculations
%                     % Here we use a simplified approximation
%                     A = (1 / r) * (1 + 0.5 * (a + b + c) / r); % Simplified
%                     B = -1 / (6 * r^2);
%                     rpy_tensor(i, j, :, :) = A * eye(3) + B * (eye(3) - ones(3, 3));
%                 end
%             end
%         end
%     end
% end


% function R = ellipsoidal_hydrodynamic_tensor(a, b, c, r)
%     % Compute the hydrodynamic resistance matrix for ellipsoidal particles
%     % a, b, c: semi-axes of the ellipsoid
%     % r: separation distance between the centers of the ellipsoids
% 
%     % Parameters
%     eta = 1; % Fluid viscosity, normalized for simplicity
% 
%     % Compute components of the resistance matrix
%     % Note: This is a simplified version; exact computations require more complex formulas
%     A = (1 / (6 * pi * eta)) * (1 / (1 + (a^2 / r^2)));
%     B = (1 / (6 * pi * eta)) * (1 / (1 + (b^2 / r^2)));
%     C = (1 / (6 * pi * eta)) * (1 / (1 + (c^2 / r^2)));
% 
%     % Resistance matrix
%     R = A * eye(3) + B * (ones(3) - eye(3)) + C * (eye(3) - ones(3));
% end

