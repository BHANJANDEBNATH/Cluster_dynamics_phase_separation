
% radius of gyration calculation 
% incorporating periodic boundaries

function [rad_gyr,unwrapped_positions] = radius_gyration_pbc(ids, positions, XX,YY, dia_disk)

N = length(ids);
unwrapped_positions = unwrap_positions_2D(ids,positions, XX, YY, dia_disk);
com = mean(unwrapped_positions, 1) ;
dis_sq = sum((unwrapped_positions - com).^2, 2);

rad_gyr = sum(dis_sq) / N;
rad_gyr = sqrt(rad_gyr);

end



