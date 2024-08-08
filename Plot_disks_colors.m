function Plot_disks_colors(No_disks, dia_disk, disks_stat, XX, YY, bondpos)

clf;
hold on;

rectangle('Position',[0 0 XX YY],'LineWidth',3)
daspect([1 1 1]);

res = 600;

for i = 1:No_disks
     
     
     if mean(disks_stat(i, bondpos(1,1):bondpos(1,length(bondpos)))) == 0
        p = nsidedpoly(1000, 'Center', [disks_stat(i,2) disks_stat(i,3)], 'Radius', dia_disk/2);
        plot(p, 'FaceColor','#D95319','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     else
         p = nsidedpoly(1000, 'Center', [disks_stat(i,2) disks_stat(i,3)], 'Radius', dia_disk/2);
         plot(p, 'FaceColor','#7E2F8E','LineWidth',1.5,'EdgeColor','black','FaceAlpha', 1)
     end
     
     %text(disks_stat(i,2), disks_stat(i,3), 1, string(disks_stat(i,1)))
    

end

    hold off;
    drawnow;



end