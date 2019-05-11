function [Image2d]=plot7_image(laVar,data,i_ok,CLim,nomJour,L,La,Dim)

minlat=min(La(:));
minlon=min(L(:));
maxlat=max(La(:));
maxlon=max(L(:));

Line = Dim(1);
Colum = Dim(2);
Cdirect=nan*ones(Line,Colum);

data(data>=CLim(2))=CLim(2);
data(data<=CLim(1))=CLim(1);
Cdirect(i_ok)=data;

Cdirect=rot90(Cdirect,1);

Image2d=Cdirect;

m_proj('miller','lat',[minlat maxlat],'lon',[minlon maxlon]);
m_coast('patch',[0.9999 0.9999 0.9999],'edgecolor','none');
m_contourf(L,La,flipud(Cdirect));
m_grid('box','fancy','tickdir','out');
colormap('jet'); 
caxis(CLim);
if strcmp(laVar(end-2:end),'SSS')
%     caxis([0 48]);
    colormap(jet(48));
elseif strcmp(laVar(end-2:end),'SST')
%     caxis([20 35]);
    colormap(jet(16));
elseif strcmp(laVar(end-2:end),'ADT')
%     caxis([-0.3 1.2]);
    colormap(jet(16));
end
colorbar;

hold on

xlabel('longitude')
ylabel('latitude')

h = title([laVar, ' du ', nomJour(end-9:end)],'fontsize',16);
set(h,'interpreter','none')


FileName=[pwd, '/figures/',laVar,'/',nomJour];
fig = gcf; fig.PaperUnits = 'inches';
fig.PaperPosition = [0 2.5 8.5 3.2];
print(FileName,'-dpng','-r0')
close
