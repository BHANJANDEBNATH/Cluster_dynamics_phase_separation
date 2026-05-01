
function unwrapped_positions = unwrap_positions_2D(ids,positions, XX, YY, dia_disk)

box_size = [XX YY];
N = length(ids); 
unwrapped_positions = positions; 

% finding out whether cluster lies on a periodic boundary
condition = [];
XL = [];
XR = [];
YB = [];
YT = [];
for i = 1:size(positions,1)
    ax = positions(i,1);
    ay = positions(i,2);
    xl = ax - 0;
    XL = [XL; xl];

    xr = XX - ax;
    XR = [XR; xr];

    yb = ay - 0;
    YB = [YB; yb];

    yt = YY - ay;
    YT = [YT; yt];
end
XL_res = XL(XL <= dia_disk); 
XR_res = XR(XR <= dia_disk);
YB_res = YB(YB <= dia_disk); 
YT_res = YT(YT <= dia_disk);

if isempty(XL_res) && isempty(XR_res) 
    condition_x = 1;
elseif ~isempty(XL_res) && isempty(XR_res)
    condition_x = 1; 
elseif isempty(XL_res) && ~isempty(XR_res)
    condition_x = 1; 
elseif ~isempty(XL_res) && ~isempty(XR_res)
    condition_x = 0;
end

if isempty(YB_res) && isempty(YT_res)
    condition_y = 1;
elseif ~isempty(YB_res) && isempty(YT_res)
    condition_y = 1; 
elseif isempty(YB_res) && ~isempty(YT_res)
    condition_y = 1; 
elseif ~isempty(YB_res) && ~isempty(YT_res)
    condition_y = 0;
end


if condition_x == 1 && condition_y == 1    % it implies that cluster does not cross periodic boundary
    unwrapped_positions = positions;
else
    % Apply unwrapping based on the minimum image convention
    for i = 1:N
        for j = 1:2 
            delta = positions(i,j) - positions(1,j);
            if delta > box_size(j)/2
                unwrapped_positions(i,j) = positions(i,j) - box_size(j); % Particle crossed right boundary
            elseif delta < - box_size(j)/2
                unwrapped_positions(i,j) = positions(i,j) + box_size(j); % Particle crossed left boundary
            else
                unwrapped_positions(i,j) = positions(i,j);               % No boundary crossing
            end
        end
    end

end


end