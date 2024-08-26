function [new_pos_s] = BoundaryCondition_disks(new_pos_s, XX, YY)



% periodic BC in X direction
if new_pos_s(1) < 0
   new_pos_s(1) = new_pos_s(1) + XX; 
end
if new_pos_s(1) > XX
   new_pos_s(1) = new_pos_s(1) - XX;
end


% periodic BC in Y direction
if new_pos_s(2) < 0
   new_pos_s(2) = new_pos_s(2) + YY; 
end
if new_pos_s(2) > YY
   new_pos_s(2) = new_pos_s(2) - YY;
end


end