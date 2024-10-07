

function [update_coord] = translate_rotate_cluster(ids,unwrapped_positions,translate,rotate,pbc_xmin,pbc_xmax,pbc_ymin,pbc_ymax);

update_coord = zeros(size(unwrapped_positions));

unwrapped_tr = zeros(size(unwrapped_positions));
unwrapped_rot = zeros(size(unwrapped_positions));

cos_theta = cos(rotate);
sin_theta = sin(rotate);

% first, perform translation operation
for i = 1:size(unwrapped_positions,1)
    unwrapped_tr(i,1) = unwrapped_positions(i,1) + translate(1);
    unwrapped_tr(i,2) = unwrapped_positions(i,2) + translate(2);
end

com = sum(unwrapped_tr, 1) / length(ids);

% second, perform rotation
for i = 1:size(unwrapped_positions,1)
    % Step 1: Translate point to center of rotation
    x_prime = unwrapped_tr(i,1) - com(1);
    y_prime = unwrapped_tr(i,2) - com(2);

    % Step 2: Rotate point using the rotation matrix
    x_new = x_prime * cos_theta - y_prime * sin_theta;
    y_new = x_prime * sin_theta + y_prime * cos_theta;
    
    % Step 3: Translate point back to original position
    unwrapped_rot(i,1) = x_new + com(1);
    unwrapped_rot(i,2) = y_new + com(2);
end

% third, apply pbc to wrap coordinates inside simulation box
for i = 1:size(unwrapped_positions,1)
    x_tr_rot = unwrapped_rot(i,1);
    y_tr_rot = unwrapped_rot(i,2);

    if x_tr_rot < pbc_xmin 
        x_adjusted = x_tr_rot + pbc_xmax;
    elseif x_tr_rot > pbc_xmax
        x_adjusted = x_tr_rot - pbc_xmax;
    else
        x_adjusted = x_tr_rot;
    end
    
    if y_tr_rot < pbc_ymin 
        y_adjusted = y_tr_rot + pbc_ymax;
    elseif y_tr_rot > pbc_ymax
        y_adjusted = y_tr_rot - pbc_ymax;
    else
        y_adjusted = y_tr_rot;
    end

    update_coord(i, :) = [x_adjusted y_adjusted];
end

           
end

