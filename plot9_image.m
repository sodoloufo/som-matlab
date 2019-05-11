function [Image2d]=plot9_image(data,i_ok,CLim,nomJour,L,La,Dim)

minlat=min(La(:));
minlon=min(L(:));
maxlat=max(La(:));
maxlon=max(L(:));

Line = Dim(1);
Colum = Dim(2);

data(data>=CLim(2))=CLim(2);
data(data<=CLim(1))=CLim(1);

Cdirect=reshape(data,Colum,Line);


%imagesc(linspace(minlon,maxlon,Colum),linspace(maxlat,minlat,Line), Cdirect,CLim);
% %

m_proj('miller','lat',[minlat maxlat],'lon',[minlon maxlon]);
m_coast('patch',[0.9999 0.9999 0.9999],'edgecolor','none');
m_contourf(L,La,flipud(Cdirect),'LineStyle','none');
%imagesc(linspace(minlon,maxlon,Colum),linspace(maxlat,minlat,Line), Cdirect,CLim);
m_grid('box','fancy','tickdir','out');
colormap('jet'); 
% cmap=hsv(4);
% colormap(cmap);
set(colorbar,'YTick',[0 1 2 3],'YTicklabel',{'0','1','2','3'})
caxis(CLim);

colorbar;
axis image;
axis xy;
hold on 

xlabel('longitude')
ylabel('latitude')
axis image

h = title(['Les Classes du ', nomJour(end-9:end)],'fontsize',18);
set(h,'interpreter','none')

FileName=[pwd, '/figures/3_Classes/',nomJour];
fig = gcf; fig.PaperUnits = 'inches';
fig.PaperPosition = [0 2.5 8.5 3.2];
print(FileName,'-dpng','-r0')
close
