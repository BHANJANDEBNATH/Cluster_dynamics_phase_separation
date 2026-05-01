
% create neighborlist
function [neighborList] = NeighborListCreation(particles, boxSize, cutoff)
    
    N = size(particles, 1);   % Number of particles
    Lx = boxSize(1);          % Box size in x direction
    Ly = boxSize(2);          % Box size in y direction
    
    % Number of cells in each direction
    cellSize = cutoff;      % Cell size should be at least equal to cutoff
    numCellsX = floor(Lx/cellSize);
    numCellsY = floor(Ly/cellSize);
    
    % Create an empty cell array for each cell
    cells = cell(numCellsX, numCellsY);
    
    % Assign each particle to a cell
    particleCells = zeros(N, 2);
    for i = 1:N
        x = particles(i, 1);
        y = particles(i, 2);
        
        cellX = min(floor(x / cellSize) + 1, numCellsX);
        cellY = min(floor(y / cellSize) + 1, numCellsY);
        
        cells{cellX, cellY} = [cells{cellX, cellY}, i];  % storing particle ids in cells
        particleCells(i, :) = [cellX, cellY];            % storing cell [rowid, columnid] corresponding to each particle
    end
    
    % Initialize neighbor list
    neighborList = cell(N, 1);
    
    % Define neighboring cell shifts to check (-1, 0, +1 in each direction)
    neighborShifts = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 0; 0, 1; 1, -1; 1, 0; 1, 1];
    
    % Build the neighbor list for each particle
    for i = 1:N
        x_i = particles(i, 1);
        y_i = particles(i, 2);
        cellX = particleCells(i, 1);
        cellY = particleCells(i, 2);
        
        % Loop over neighboring cells
        for shiftIdx = 1:size(neighborShifts, 1)
            neighborCellX = mod(cellX + neighborShifts(shiftIdx, 1) - 1, numCellsX) + 1;
            neighborCellY = mod(cellY + neighborShifts(shiftIdx, 2) - 1, numCellsY) + 1;
            
            % Check particles in the neighboring cell
            neighbors = cells{neighborCellX, neighborCellY};
            for j = neighbors
                if i ~= j
                    
                    % dumping only those disk ids in the neighboring cells inclding periodic images
                    % who are  inside cutoff radius 
                    % % Compute distance considering periodic boundaries
                    dx = particles(j, 1) - x_i;
                    dy = particles(j, 2) - y_i;
                    dx = dx - round(dx / Lx) * Lx;  % Periodic boundary condition in x
                    dy = dy - round(dy / Ly) * Ly;  % Periodic boundary condition in y

                    % Calculate distance
                    dist = sqrt(dx^2 + dy^2);

                    % If within the cutoff distance, add to the neighbor list
                    if dist <= cutoff
                        neighborList{i} = [neighborList{i}, j];
                    end
                    
                end
            end
        end
    end


end
