

function Plot_disks_colors_with_images(No_disks, dia_disk, disks_coord_plot, XX, YY, bondpos)

clf;
hold on;

rectangle('Position',[(0 - XX) (0) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0 - XX) (0 - YY) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0 - XX) (0 + YY) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0 + XX) (0) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0 + XX) (0 - YY) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0 + XX) (0 + YY) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0) (0 - YY) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0) (0 + YY) (XX) (YY)],'LineWidth',3)
daspect([1 1 1]);

rectangle('Position',[(0) (0) (XX) (YY)],'LineWidth',3, 'EdgeColor','red')
daspect([1 1 1]);

res = 600;

for i = 1:No_disks
     
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [disks_coord_plot(i,2) disks_coord_plot(i,3)], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [disks_coord_plot(i,2) disks_coord_plot(i,3)], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end
     
     %text(disks_coord_plot(i,2), disks_coord_plot(i,3), 1, string(disks_coord_plot(i,1)))
    
end

for i = 1:No_disks
    x2 = disks_coord_plot(i,2) + XX; 
    y2 = disks_coord_plot(i,3);
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x2 y2], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x2 y2], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end
end

for i = 1:No_disks
    x3 = disks_coord_plot(i,2) + XX; 
    y3 = disks_coord_plot(i,3) + YY;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x3 y3], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x3 y3], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end 
end

for i = 1:No_disks
    x4 = disks_coord_plot(i,2) + XX; 
    y4 = disks_coord_plot(i,3) - YY;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x4 y4], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x4 y4], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end    
end

for i = 1:No_disks
    x5 = disks_coord_plot(i,2) ; 
    y5 = disks_coord_plot(i,3) + YY;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x5 y5], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x5 y5], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end    
end

for i = 1:No_disks
    x6 = disks_coord_plot(i,2) ; 
    y6 = disks_coord_plot(i,3) - YY;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x6 y6], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x6 y6], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end      
end

for i = 1:No_disks
    x7 = disks_coord_plot(i,2) - XX ; 
    y7 = disks_coord_plot(i,3) ;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x7 y7], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x7 y7], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end     
end

for i = 1:No_disks
    x8 = disks_coord_plot(i,2) - XX ; 
    y8 = disks_coord_plot(i,3) + YY;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x8 y8], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x8 y8], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end     
end

for i = 1:No_disks
    x9 = disks_coord_plot(i,2) - XX ; 
    y9 = disks_coord_plot(i,3) - YY;
     if mean(disks_coord_plot(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [x9 y9], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [x9 y9], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#77AC30','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
         %plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end    
end

hold off;
drawnow;

end