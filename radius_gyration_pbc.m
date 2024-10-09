
% radius of gyration calculation 
% incorporating periodic boundaries

function [rad_gyr,unwrapped_positions] = radius_gyration_pbc(ids, positions, XX,YY, dia_disk)

N = length(ids);

unwrapped_positions = unwrap_positions_2D(ids,positions, XX, YY, dia_disk);

com = mean(unwrapped_positions, 1) ;

rad_gyr = 0;

for i = 1:N
    % Compute the shortest distance between the particle and COM under PBC
    delta = unwrapped_positions(i,:) - com;
    
    rad_gyr = rad_gyr + sum(delta.^2);
end
    
rad_gyr = rad_gyr / N;
rad_gyr = sqrt(rad_gyr);

end



