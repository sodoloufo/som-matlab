clear all ; close all; clc ;
addpath ('som');
addpath ('som\dijkstradir');
cnames = {'sss','sst','adt'};

load results_app2_40_30.mat 
load  Tot_SSS.mat; 
load Tot_SST.mat;                         
                                              
load Tot_ADT.mat;                            

load Tot_EWSS.mat;                             
 
load Tot_NSSS.mat;                             
for i=1:length(Tot_EWSS)
norm_vent(i)= norm([Tot_EWSS(i),Tot_NSSS(i)]);
end
normmm=norm_vent';
base_cr = [Tot_SSS Tot_SST Tot_ADT];

%--------------------------------------------------------------
% M?lange aleatoire des donnees
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
data = base_test_mm(Appstart:Appstart+Nens, :);  

sDapp = som_data_struct(data,'comp_names',cnames);
essai=sMap.topol.msize; nbneuronS=essai(1)*essai(2); 
nlmap  = essai(1); ncmap = essai(2) ; % Nombre de lignes et nombre de colonnes 

nb_ref = nlmap * ncmap;  % Nombre de referents : nbre de lignes x nbre de colonnes
msize  = [nlmap ncmap];  % Dimension de la carte
%% 2) Classification hierarchique des r?f?rents %%

carte_ref_et_donn(sMap,data) % Deploiement des r?f?rents (de la carte) sur les donn?es.

class_data = classification(sMap,data);
class_dat=class_data; % classe de donnees
bmus_1=bmus_smap3(find(class_dat==1)); % Pour trouver tous les bmus correspondant ? la classe 1
val_b1=[];
for i=1:length(bmus_1)                 % Pour trouver les neurones de la carte smap3 correspondant ? bmus de la classe 1
    val_b=smap3.codebook(bmus_1(i),:)
    val_b1=[val_b1;val_b];
end
mean_b1=mean(val_b1);                  % moyenne de la sst et de la chloro