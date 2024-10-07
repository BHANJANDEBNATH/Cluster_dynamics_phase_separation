
% radius of gyration calculation 
% incorporating periodic boundaries

function [rad_gyr,unwrapped_positions] = radius_gyration_pbc(ids, positions, XX,YY)

N = length(ids);
box_size = [XX YY];

unwrapped_positions = unwrap_positions_2D(ids,positions, XX, YY);
com = sum(unwrapped_positions, 1) / N;

rad_gyr = 0;

for i = 1:N
    % Compute the shortest distance between the particle and COM under PBC
    delta = unwrapped_positions(i,:) - com;
    delta = delta - box_size .* round(delta ./ box_size);  % Apply periodic boundary condition
    
    rad_gyr = rad_gyr + sum(delta.^2);
end
    
rad_gyr = rad_gyr / N;
rad_gyr = sqrt(rad_gyr);

end

function unwrapped_positions = unwrap_positions_2D(ids,positions, XX, YY)

box_size = [XX YY];
N = length(ids); 
unwrapped_positions = positions; 
    
    % Apply unwrapping based on the minimum image convention
    for i = 2:N
        for j = 1:2 
            delta = positions(i,j) - positions(i-1,j);
            if delta > box_size(j)/2
                unwrapped_positions(i,j) = positions(i,j) - box_size(j); % Particle crossed right boundary
            elseif delta < -box_size(j)/2
                unwrapped_positions(i,j) = positions(i,j) + box_size(j); % Particle crossed left boundary
            else
                unwrapped_positions(i,j) = positions(i,j); % No boundary crossing
            end
        end
    end


end


