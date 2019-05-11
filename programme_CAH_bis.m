clear all ; close all; clc ;
addpath ('som');
addpath ('som\dijkstradir');
addpath('base_traite');

cnames = {'SSS', 'SST','ADT'};

load results_app2.mat
load Tot_SSS.mat;
load Tot_SST.mat;
load Tot_ADT.mat;

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


nlmap  = sMap.topol.msize(1); ncmap = sMap.topol.msize(2) ; % Nombre de lignes et nombre de colonnes
nb_ref = nlmap * ncmap;  % Nombre de referents : nbre de lignes x nbre de colonnes
msize  = sMap.topol.msize;  % Dimension de la carte
nomCarte = [num2str(msize(1)),'vx',num2str(msize(2))];
%% 2) Classification hierarchique des r?f?rents %%

%% matrice de distance euclidienne entre les r?ferents 2 ? 2 %%
% rng('default');  % For reproducibility
% eva = evalclusters(sMap.codebook,'linkage','CalinskiHarabasz','KList',[1:6])
%
D = pdist(sMap.codebook, 'euclidean');

% db_index
% [t,r] = db_index(sMap.codebook,600)
%
%%%%%%%% Classification %%%%%

Z = linkage(D,'complete');

figure1 = figure('color',[1 1 1]);
bar (Z(nlmap*ncmap-1:-1:nlmap*ncmap-20,3));
title('Diagramme en bar')

%%% Affichage de l?arbre hi?rarchique des clusters de donnees %%%

figure2= figure('color',[1 1 1]);
dendrogram(Z,nlmap*ncmap);
hold on;
plot ([0,300],[1.5,1.5]);
%dendrogram(Z,'ColorThreshold','default');
title('Arbre hierarchique des clusters de donn?es')
%%%repr?sentation des 20 derni?res distances d?agr?gation des clusters
%%%(3?me colonne de Z)%%



%%%% les classes des referents %%%
nb_class = 3
class_ref = cluster (Z,nb_class);
save class_ref class_ref

%%% Presentation de la classification obtenue sur la carte topologique%%
figure3= figure('color',[1 1 1]);
LL=[1 2 3]; sMap = som_label(sMap, 'clear','all');
sMap = som_label(sMap,'add',[1:nb_ref]',int2str([1:nb_ref]'));
sMap = som_autolabel(sMap, sDapp,'freq');
som_cplane(sMap, class_ref,1);

colormap('jet');
set(colorbar,'YTick',[0 1 2 3],'YTicklabel',{'0','1','2','3'})
title('Classification hierarchique sur les referents','fontsize',14)
FileName=[pwd, '/figures/cah/CAH_Jet_3', nomCarte];
print(FileName,'-dpng','-r0')

%% %%  Affectation de chaque BMUS a une classe donnee

indneurcarte=find(sMap.codebook(:,1));
indclass1=find(class_ref==1);
indclass2=find(class_ref==2);
indclass3=find(class_ref==3);
% indclass4=find(class_ref==4)
% indclass5=find(class_ref==5)

neur_class1=indneurcarte(indclass1);
neur_class2=indneurcarte(indclass2);
neur_class3=indneurcarte(indclass3);
% neur_class4=indneurcarte(indclass4);
% neur_class5=indneurcarte(indclass5);

weiclasse1 = sMap.codebook(indclass1,:);
weiclasse2 = sMap.codebook(indclass2,:);
weiclasse3 = sMap.codebook(indclass3,:);
% weiclasse4 = sMap.codebook(indclass4,:);
% weiclasse5 = sMap.codebook(indclass5,:);



save smap_class weiclasse1 weiclasse2 weiclasse3  neur_class1 neur_class2 neur_class3 
% weiclasse4 neur_class4

bmus_1=bmus(find(class_ref==1)); % Pour trouver tous les bmus correspondant ? la classe 1
val_b1=[];
for i=1:length(bmus_1)
    val_b=find(bmus==indclass1(i))
    val_b1=[val_b1;val_b];
end
bmus_1=val_b1;

bmus_2=bmus(find(class_ref==2)); % Pour trouver tous les bmus correspondant ? la classe 1
val_b2=[];
for i=1:length(bmus_2)
    val_b=find(bmus==indclass2(i))
    val_b2=[val_b2;val_b];
end
bmus_2=val_b2;

bmus_3=bmus(find(class_ref==3)); % Pour trouver tous les bmus correspondant ? la classe 1
val_b3=[];
for i=1:length(bmus_3)
    val_b=find(bmus==indclass3(i))
    val_b3=[val_b3;val_b];
end
bmus_3=val_b3;

% bmus_4=bmus(find(class_ref==4)); % Pour trouver tous les bmus correspondant ? la classe 1
% val_b4=[];
% for i=1:length(bmus_4)
%     val_b=find(bmus==indclass4(i))
%     val_b4=[val_b4;val_b];
% end
% bmus_4=val_b4;
%
% bmus_5=bmus(find(class_ref==5)); % Pour trouver tous les bmus correspondant ? la classe 1
% val_b5=[];
% for i=1:length(bmus_5)
%     val_b=find(bmus==indclass5(i))
%     val_b5=[val_b5;val_b];
% end
% bmus_5=val_b5;

% figure;
% imagesc(reshape(class_ref,nlmap,ncmap))

%figure;
%%observation classe d appartenance
classe_dap = class_ref(bmus)

for i=1:3
    eval(['indapp' num2str(i) ' = find(classe_dap==' num2str(i) ');']);
end

% for i= 1:5
%  for k = 1:size(eval(['indapp' num2str(i) ]))
%         eval(['class_dap' num2str(i) ' = sDapp.data(indapp' num2str(i) ',' k ');']);
%  end
% end
%
%pour retrouver les classess
class_dap1 = sDapp.data(indapp1,:)
class_dap2 = sDapp.data(indapp2,:)
class_dap3 = sDapp.data(indapp3,:)
% class_dap4 = sDapp.data(indapp4,:)

save class_ref class_dap1 class_dap2 class_dap3 class_dap3
% class_dap4
%  base = ["indapp1" ,"indapp2", "indapp3" ,"indapp4", "indapp5"]


%%Avec boucle
for i=1:nlmap * ncmap
    eval(['Captes' num2str(i) ' = find(bmus ==' num2str(i) ');']);
end

%%
% figure2 = figure('color',[1 1 1]);hold on
% imagesc(classe1(:,1),25,36)
% imagesc(classe2(:,2),25,36)
%
% data1 = [mean(class_dap1(:,1)),mean(class_dap1(:,2)),mean(class_dap1(:,3))];
% figure ;
% c = categorical(cnames);
% bar(data1,0.4);
% ylabel('Moyennes des Valeurs des variables');
% xlabel('sss,sst,adt');
% title('Presentation de l histogramme des moyennes classe 1')
%
%
% data2 = [mean(class_dap2(:,1)),mean(class_dap2(:,2)),mean(class_dap2(:,3))];
% figure ;
% c = categorical(cnames);
% bar(c,data2,0.4);
% ylabel('Moyennes des Valeurs des variables');
% xlabel('sss,sst,adt');
% title('Presentation de l histogramme des moyennes classe 2')
%
%
% data3 = [mean(class_dap3(:,1)),mean(class_dap3(:,2)),mean(class_dap3(:,3))];
%
% figure ;
% c = categorical(cnames);
% bar(c,data3,0.4);
% ylabel('Moyennes des Valeurs des variables');
% xlabel('sss,sst,adt');
% title('Presentation de l histogramme des moyennes classe 3')


% data4 = [sum(class_dap4(:,1)),sum(class_dap4(:,2)),sum(class_dap4(:,3))];
%
% figure ;
% c = categorical(cnames);
% bar(c,data4,0.4);
% ylabel('pixel');
%%
figure3 = figure('color',[1 1 1]);hold on
plot(class_dap1(:,1),class_dap1(:,2),'c+');
plot(class_dap2(:,1),class_dap2(:,2),'y*');
plot(class_dap3(:,1),class_dap3(:,2),'ro');
%plot(class_dap4(:,1),class_dap4(:,2),'ko')


for i=1:nb_ref
    exemples_captes = find(bmus== i);%indice de la donnee captee par le neurone x
    exemples_affecte = sMap.codebook(i,2);%affecte a chaque donnee son poids correspondant sur smap
end