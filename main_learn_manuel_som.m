clear all ; close all; clc ;

addpath ('som');
addpath ('som\dijkstradir');
                           
load Masque.mat 
load nys.mat 
load nxs.mat
load xlos.mat
load xlas.mat

load('Tot_SSS.mat') 
load('Tot_SST.mat')                          
load('Tot_ADT.mat') 

base_cr = [Tot_SSS Tot_SST Tot_ADT];

%--------------------------------------------------------------
% M??lange aleatoire des donnees
iperm    = randperm(size(base_cr,1));
base_mix = base_cr(iperm,:);
%save base_mix base_mix % si on veut
%
N=size(base_mix,1);

for i=1:3
[base_test_mm(:,i),  mini(:,i), maxi(:,i) ]= normal_min_max(base_mix(:,i));
end
save mini mini 
save maxi maxi 
clear base_mix  Tot_* base_cr
%
%--------------------------------------------------------------
% Constitution des ensembles (App, Val et Test)
%
% Positionnement d'un indice de debut de selection pour chaque ensemble
Appstart = 1;
Valstart = floor(N/3);  
Teststart= floor(N*2/3);
Nens     = ceil(N*0.02);    % Taille des ensembles : 2%
%
% Ensemble d'apprentissage
Xa = base_test_mm(Appstart:Appstart+Nens, :);  

% Ensemble de validation
Xv = base_test_mm(Valstart:Valstart+Nens, :);   

% Ensemble de Test
Xt = base_test_mm(Teststart:Teststart+Nens, :); 

% structure de donnees
% --------------------
% cnames = {'sss', 'sst','adt'};
% structData= som_data_struct(Xa,'comp_names',cnames);
% sData=som_normalize(structData,'range');

cnames = {'SSS', 'SST','ADT'};
% structData= som_data_struct(Xa,'comp_names',cnames);
% sData=som_normalize(structData,'range');
sData= som_data_struct(Xa,'comp_names',cnames);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection ACP
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Pd,V,me]=pcaproj(sData,2);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation de la carte
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msize =[33 30];
nomCarte = [num2str(msize(1)),'vx',num2str(msize(2))];
insize  = size(sData.data,2);
lattice = 'hexa';     %% 'hexa' ou 'rect'
shape   = 'sheet';    %% 'sheet', 'cyl', ou 'toroid'
sMap = som_map_struct(insize,'msize',msize, lattice, shape);
% Affichage des donn??es ?? apprendre
sMap_initiale = som_lininit(sData, sMap);

figure
plot_data_proj(sMap_initiale,Pd,V,me);
title('Etat initial','fontsize',16);
FileName=[pwd, '/figures/apprentissage/Etat_initial',nomCarte];
print(FileName,'-dpng','-r0')
close

sMap = sMap_initiale;
figure
%% Premier apprentissage : 
%%% --> Variation rapide de la temperature de 6 a 2 en 90 iterations
epochs     = 60; % 120 %150
radius_ini =9; %10;
%radius_fin = max(1,radius_ini/4);
radius_fin = 1;
Neigh      = 'gaussian'; %% 'gaussian', 'cutgauss', 'bubble' ou 'ep'
tr_lev     = 3;
D=sData(1:10:end,:);
[sMap,sTruct] = som_batchtrain(sMap,D,'trainlen',epochs,...
     'radius_ini',radius_ini,'radius_fin',radius_fin, ...
	 'neigh',Neigh,'tracking',tr_lev);
% Visualisatin des donnees et de la Map
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMap1 = sMap;
som_show(sMap);
print -djpeg sMap1erapp_[33 30]
figure
plot_data_proj(sMap,Pd,V,me);
title('Etat du premier apprentissage');
print -djpeg Etat_premier_apprentissage_[33 30]
%break;
figure;
% affichage de la carte
% --------------------
%som_show(sMap) 
som_show(sMap,'umat','all')
%U=som_umat(sMap);
title('U-matrix');
print -djpeg umat2_[33 30]
figure;
som_cplane('hexa',[33 30],sMap.codebook(:,1));
colorbar;
title('sss');
print -djpeg sss1_[33 30]
figure;
som_cplane('hexa',[33 30],sMap.codebook(:,2));
colorbar;
title('sst');
print -djpeg sst1_[33 30]
figure;
som_cplane('hexa',[33 30],sMap.codebook(:,3));
colorbar;
title('\adt');
print -djpeg adt1_[33 30]
% figure;
% som_cplane('hexa',[33 30],sMap.codebook(:,4));
% colorbar;
% title('\wind');
% print -djpeg adt1_[33 30]

%%%%calcul des erreurs de quantification et topographique  
[mqe,tge] = som_quality(sMap,D);
fprintf(1,'Final quantization error: %5.3f\n',mqe)
fprintf(1,'Final topographic error:  %5.3f\n',tge)

%%%%%%%%%%%%% WFEP %%%%%%%%%%%%%%%%%%%%%%
bmus = som_bmus(sMap,D);
hits = som_hits(sMap,D);

save results_app1_33_30 bmus hits sMap;
%break;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% Deuxieme apprentissage :
%%% --> Variation fine de la temperature de 1 a 0.7 en 200 iterations
epochs     = 90; %300  %200
radius_ini = 1;  % 2
radius_fin = 0.1;   %0.1
D=sData(1:10:end,:);
[sMap,sTruct] = som_batchtrain(sMap,D,'trainlen',epochs,...
     'radius_ini',radius_ini,'radius_fin',radius_fin, ...
	 'neigh',Neigh,'tracking',tr_lev);

%save workspace;
% Visualisatin des donnees et de la Map
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sMap2 = sMap;
figure;
som_show(sMap);
FileName=[pwd, '/figures/apprentissage/sMap2app', nomCarte];
print(FileName,'-dpng','-r0')
close

figure
plot_data_proj(sMap,Pd,V,me);
title('Etat du deuxieme apprentissage','fontsize',16);
FileName=[pwd, '/figures/apprentissage/Etat_deuxieme_apprentissage', nomCarte];
print(FileName,'-dpng','-r0')


% affichage de la carte
% --------------------
figure;
% affichage de la carte
som_show(sMap,'umat','all')
title('U-matrix');
FileName=[pwd, '/figures/apprentissage/umat2', nomCarte];
print(FileName,'-dpng','-r0')
close

figure;
som_cplane('hexa',sMap.topol.msize,sMap.codebook(:,1));
colorbar;
title('SSS');
FileName=[pwd, '/figures/apprentissage/sss2', nomCarte];
print(FileName,'-dpng','-r0')
close

figure;
som_cplane('hexa',sMap.topol.msize,sMap.codebook(:,2));
colorbar;
title('SST');
FileName=[pwd, '/figures/apprentissage/sst2', nomCarte];
print(FileName,'-dpng','-r0')
close

figure;
som_cplane('hexa',sMap.topol.msize,sMap.codebook(:,3));
colorbar;
title('ADT');
FileName=[pwd, '/figures/apprentissage/adt2', nomCarte];
print(FileName,'-dpng','-r0')
close

%%%%calcul des erreurs de quantification et topographique  
[mqe,tge] = som_quality(sMap,D);
fprintf(1,'Final quantization error: %5.3f\n',mqe)
fprintf(1,'Final topographic error:  %5.3f\n',tge)

%%%%%%%%%%%%% WFEP %%%%%%%%%%%%%%%%%%%%%%
bmus = som_bmus(sMap,D);
hits = som_hits(sMap,D);

save results_app2 bmus hits sMap;