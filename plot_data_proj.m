function plot_data_proj(sM,Pd,V,me,labscol)
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fonction plot_data_proj
%%%
%%% [Pd,V,me]=pcaproj(sD,2);
%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_data = size(Pd,1);
% affiche les donnees (la projection)
if nargin < 5
  h=som_grid('rect',[n_data 1], 'markersize',3,'Line','none','Coord',Pd);
else
  h=som_grid('rect',[n_data 1], 'markersize',3,'markercolor',labscol, ...
	     'Line','none','Coord',Pd);
end
% affiche aussi la Map avec la projection ACP des Donnees
hold on, grid on, box on
h=som_grid(sM,'Coord',pcaproj(sM,V,me),'marker','none', ...
	   'linecolor',[0.5 0.5 0.5]);
hold off, axis on
%axis image
title('ACP Projection');

return
