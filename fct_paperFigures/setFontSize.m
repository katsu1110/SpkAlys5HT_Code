function setFontSize( ax, fsz)
%sets the font size to fsz and the label font size to fsz + x
global fsz_add;

ax.FontSize = fsz;
ax.XLabel.FontSize = fsz+fsz_add;
ax.YLabel.FontSize = fsz+fsz_add;


ax.Box = 'off';  
ax.TickDir = 'out';           

end

