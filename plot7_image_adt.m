function [Image2d]=plot7_image_adt(numClass,data,i_ok,CLim,nomJour,L,La,Dim)

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

Image2d=Cdirect;
Image2d=flipud(Image2d');
ii=find(Image2d<0.002);
Image2d(ii)=NaN;
iii=find(Image2d>=0.002);
Image2d(iii)=numClass;
Image2d=Image2d(:);
Image2d=reshape(Image2d,Colum,Line);
% Cdirect=rot90(Cdirect,1);
% 
% 
% m_proj('miller','lat',[minlat maxlat],'lon',[minlon maxlon]);
% m_coast('patch',[0.9999 0.9999 0.9999],'edgecolor','none');
% m_contourf(L,La,flipud(Cdirect));
% m_grid('box','fancy','tickdir','out');
% colormap('jet'); colorbar;
% caxis(CLim);
% 
% hold on
% 
% xlabel('longitude')
% ylabel('latitude')
% h = title(['ADT Classe ',num2str(numClass),' du ', nomJour(end-9:end)],'fontsize',12);
% set(h,'interpreter','none')
% 
% 
% FileName=[pwd, '/figures/classe_var/',nomJour,'/ADT_Classe',num2str(numClass)];
% fig = gcf; fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 2.5 8.5 3.2];
% print(FileName,'-dpng','-r0')
% close
