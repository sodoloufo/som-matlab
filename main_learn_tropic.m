clear all ; close all; close all ;

addpath ('E:\mes projets\kaly_Ibra\cmos2012_2014\SOM-Toolbox_v2.1\som');
addpath ('E:\mes projets\kaly_Ibra\cmos2012_2014\SOM-Toolbox_v2.1\som\dijkstradir');
load  Tot_SSS.mat;                            
load Masque.mat 
load nys.mat 
load nxs.mat
load xlos.mat
load xlas.mat

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
% Mélange aleatoire des donnees
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

%************************************************************************************************
%                                        LA CARTE TOPOLOGIQUE
% On fait cette fois des cartes sur chaque variable SST ADT SSS
%************************************************************************************************
%
% carte topologique temporelle SST
sMapT=som_make(Xa); 
sData=som_data_struct(Xa);


Don=som_normalize(sData,'var');

%sMapT=som_make(Don);

sMapT=som_make(Don,'munits',800,'training','long','lattice','hexa');

sMapT=som_batchtrain(sMapT,Don,'radius',[10 5],'trainlen',20,'tracking',3);
sMapT=som_batchtrain(sMapT,Don,'radius',[5 3],'trainlen',50,'tracking',3);
sMapT=som_batchtrain(sMapT,Don,'radius',[3 1],'trainlen',80,'tracking',3);
sMapT=som_batchtrain(sMapT,Don,'radius',[1 0.01],'trainlen',160,'tracking',3);

%som_show(sMapT,'comp',[1:4],'umat','all');
som_show(sMapT,'umat','all');

% map size [15, 8]
% Final quantization error: 61.626
% Final topographic error: 0.024
% on utilise donc  15*8 neurones
save sMapT sMapT


nombrespace=size(Xa,2);

UT=sMapT.codebook;
bmusT=som_bmus(sMapT,Xa); % on donne a chaque pixel de la carte, la valeur de son référent
save bmusT bmusT 
for i=1:nombrespace
    XcarteT(:,i)=UT(bmusT,i);
end
essai=sMapT.topol.msize; nbneuronT=essai(1)*essai(2); clear essai


%load('Cartes_separees.mat')

%%

% % denormalisation
% %********************************
% XTd= decentred(XcarteT,moyT,ecT,1);
% clear XcarteT;
% XcarteT=XTd; clear XTd
% 
% % decentrage et dereduction de SSS     
% %*********************************
% XSd = decentred(XcarteS,moyS,ecS,1);
% clear XcarteS;
% XcarteS=XSd; clear XSd
% 
% % decentrage et dereduction de ADT  
% %*********************************
% XAd = decentred(XcarteA,moyA,ecA,1);
% clear XcarteA;
% XcarteA=XAd; clear XAd

XTd= denorm_min_max(XcarteT(:,1), mini(:,1), maxi(:,1));
clear XcarteT;
XcarteT=XTd; clear XTd


%************************************************************************************************
%                                        CLASSIFICATION
%************************************************************************************************
% SST
%*****
ZT = linkage(pdist(UT,'euclidean'),'ward'); 
% linkage cree un arbre de classification avec euclidean=norme euclidienne, ward=distance interne carree
% pdist calcule la distance euclidienne entre 2 obs

% affichage du dendrogramme, choix des classes
subplot(2,1,1)
dendrogram(ZT,nbneuronT); set(gca,'xtick',[]);% trace le dendrogramme avec nbneuronT branches sans afficher les labels ox
% Remarque: sans mettre nbneuronT, le dendrogramme elimine les branches inferieures pour clarifier la figure
title('Dendrogram Eddy','fontsize',14); 
subplot(2,1,2)
bar(ZT(nbneuronT-1:-1:nbneuronT-20,3)); axis([0 21 0 1500]);%  affiche le nb de données ds chaque cluster juste pour les 20
title('data number for the last 20 clusters','fontsize',14);
%print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_dendroT_v2018.png' 
%close
 
% on essaie 3 classes
n=3;
class_ref_T =cluster(ZT,n); % classification de ZT en n classes
som_cplane (sMapT, class_ref_T,1); %visualise la carte de Kohonen des 3 classes en T
%print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_mapT_v2018.png';
%close
% som_bmus va permettre d'associer chaque individu à son référent puis sa classe
bmusT=som_bmus(sMapT,Xa); classT=class_ref_T(bmusT);


% representation temporelle des classes pour la SST
for i=1:length(classT)
    if classT(i)==1
           classT1(i)=classT(i); classT2(i)=NaN; classT3(i)=NaN;
    elseif classT(i)==2
           classT2(i)=classT(i); classT3(i)=NaN; classT1(i)=NaN;
    else
           classT3(i)=classT(i); classT2(i)=NaN; classT1(i)=NaN;
    end
end

% Trace l occurence des classes de T en fonction du temps
%plot(dates,classT1,'c*',dates,classT2,'r*',dates,classT3,'b*');
%axis([dates(1) dates(nts) 0 n+1]);
%set(gca,'YTick',[0 1 2 3 4]); set(gca,'XTick',[dates(1):365:dates(nts)]); 
%datetick('x','yyyy','keeplimits');
%xlabel('Time','fontsize',9); ylabel('Class number','fontsize',9);
%title('Time serie of SST eddy classes','fontsize',11);
%print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_occurenceT_v2018.png';
%close
%%
% reconstruction de la carte T pour la classe 1
j1=find(classT1==1);
nbT1=length(j1);
T1=XcarteT(j1,:);

ipt=0;
for j=1:nys % car SST et SSS ont les memes dimensions en lat lon
for i=1:nxs
if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ct1(i,j,1:nbT1) = T1(ipt); % si valeur ocean du masque donc 1 on place la valeur de la classe
%if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ct1(i,j,1:nbT1) = T1(1:nbT1, ipt);
else
Ct1(i,j,1:nbT1) = NaN;
end %if
end %for ix
end %for iy
MCt1=nanmean(Ct1,3); % moyenne de la classe 1 en T

% reconstruction de la carte T pour la classe 2
j2=find(classT2==2);
nbT2=length(j2);
T2=XcarteT(j2,:);

ipt=0;
for j=1:nys
for i=1:nxs
if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ct2(i,j,1:nbT2) = T2(ipt); % si valeur ocean du masque donc 1 on place la valeur de la classe
else
Ct2(i,j,1:nbT2) = NaN;
end %if
end %for ix
end %for iy
MCt2=nanmean(Ct2,3);

% reconstruction de la carte T pour la classe 3
j3=find(classT3==3);
nbT3=length(j3);
T3=XcarteT(j3,:);

ipt=0;
for j=1:nys
for i=1:nxs
if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ct3(i,j,1:nbT3) = T3(ipt); % si valeur ocean du masque donc 1 on place la valeur de la classe
else
Ct3(i,j,1:nbT3) = NaN;
end %if
end %for ix
end %for iy
MCt3=nanmean(Ct3,3);

clear j1 j2 j3 Masque ;

%******************************************************
% affichage des moyennes des 3 classes de la SST
set (gca,'fontsize',8)
subplot(3,2,[1 2])
contourf(xlos,xlas,MCt1');
axis equal; 
axis([-70 -30 0 15])
colorbar
%xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
title('SST eddy mean class 1','fontweight','bold','fontsize',10)

subplot(3,2,[3 4])
contourf(xlos,xlas,MCt2')
axis equal; 
axis([-70 -30 0 15])
colorbar
%xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
title('SST eddy mean class 2','fontweight','bold','fontsize',10)

subplot(3,2,[5 6])
contourf(xlos,xlas,MCt3')
axis equal; 
axis([-70 -30 0 15])
colorbar
xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
title('SST eddy mean class 3','fontweight','bold','fontsize',10)

%%print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_meanclassT_v2018.png';
%%close

% %%
% %********************
% % classification ADT                
% %********************
% ZA = linkage(pdist(UA,'euclidean'),'ward');
% subplot(2,1,1)
% dendrogram(ZA,nbneuronA); set(gca,'xtick',[]);
% title('Dendrogram','fontsize',14); 
% subplot(2,1,2)
% bar(ZA(nbneuronA-1:-1:nbneuronA-20,3)); 
% title('data number for the last 20 clusters','fontsize',14);
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_dendroA_v2018.png' 
% %close
% 
% % on essaie 3 classes
% n=3;
% class_ref_A =cluster(ZA,n);% classification de ZA en n classes
% som_cplane (sMapA, class_ref_A,1); %visualise la carte de Kohonen des 3 classes en A
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_mapA_v2018.png';
% %close
% % som_bmus va permettre d'associer chaque individu à son référent et sa classe
% bmusA=som_bmus(sMapA,XA); classA=class_ref_A(bmusA);
% 
% % representation temporelle des classes pour A
% for i=1:length(classA)
%     if classA(i)==1
%            classA1(i)=classA(i); classA2(i)=NaN; classA3(i)=NaN;
%     elseif classA(i)==2
%            classA2(i)=classA(i); classA3(i)=NaN; classA1(i)=NaN;
%     else   classA(i)==3
%            classA3(i)=classA(i); classA2(i)=NaN; classA1(i)=NaN;
%     end
% end
% 
% % Trace l occurence des classes de A en fonction du temps
% plot(dates,classA1,'c*',dates,classA2,'r*',dates,classA3,'b*')
% axis([dates(1) dates(nts) 0 n+1]);
% set(gca,'YTick',[0 1 2 3 4]); set(gca,'XTick',[dates(1):365:dates(nts)]); 
% datetick('x','yyyy','keeplimits');
% title('Time serie of ADT classes','fontsize',11);
% xlabel('Time','fontsize',9); ylabel('Class numer','fontsize',9);
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_occurenceA_v2018.png';
% %close
% 
% %%
% % reconstruction de la carte A pour la classe 1
% j1=find(classA1==1);
% nbA1=length(j1);
% A1=XcarteA(j1,:);
% 
% ipt=0;
% for j=1:nys
% for i=1:nxs
% if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ca1(i,j,1:nbA1) = A1(1:nbA1, ipt); % si valeur ocean du masque donc 1 on place la valeur de la classe
% else
% Ca1(i,j,1:nbA1) = NaN;
% end %if
% end %for ix
% end %for iy
% MCa1=nanmean(Ca1,3);
% 
% % reconstruction de la carte A pour la classe 2
% j2=find(classA2==2);
% nbA2=length(j2);
% A2=XcarteA(j2,:);
% 
% ipt=0;
% for j=1:nys
% for i=1:nxs
% if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ca2(i,j,1:nbA2) = A2(1:nbA2, ipt); % si valeur ocean du masque donc 1 on place la valeur de la classe
% else
% Ca2(i,j,1:nbA2) = NaN;
% end %if
% end %for ix
% end %for iy
% MCa2=nanmean(Ca2,3);
% 
% % reconstruction de la carte A pour la classe 2
% j3=find(classA3==3);
% nbA3=length(j3);
% A3=XcarteA(j3,:);
% 
% ipt=0;
% for j=1:nys
% for i=1:nxs
% if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Ca3(i,j,1:nbA3) = A3(1:nbA3, ipt); % si valeur ocean du masque donc 1 on place la valeur de la classe
% else
% Ca3(i,j,1:nbA3) = NaN;
% end %if
% end %for ix
% end %for iy
% MCa3=nanmean(Ca3,3);
% 
% clear j1 j2 j3
% 
% % affichage des moyennes des 3 classes de la ADT
% set (gca,'fontsize',8)
% subplot(3,2,[1 2])
% contourf(xlos,xlas,MCa1',-0.03:0.005:0.03)
% axis equal; 
% axis([-70 -30 0 15]);
% colorbar
% %xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title(['ADT mean class 1'],'fontweight','bold','fontsize',10)
% 
% subplot(3,2,[3 4])
% contourf(xlos,xlas,MCa2',-0.03:0.005:0.03)
% axis equal; 
% axis([-70 -30 0 15]);
% colorbar;
% %xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title(['ADT mean class 2'],'fontweight','bold','fontsize',10)
% 
% subplot(3,2,[5 6])
% contourf(xlos,xlas,MCa3',-0.03:0.005:0.03)
% axis equal; 
% axis([-70 -30 0 15]);
% colorbar
% xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title(['ADT mean class 3'],'fontweight','bold','fontsize',10)
% 
% %print -dpng '/usr/lodyc/tropic/sa/DATA/SMOS/SMOS2016/SVD/figures/KOHONEN_V2018/Kohonen_time_eddy_meanclassA_v2018.png';
% %close
% 
% %%
% %********************
% % classification SSS                
% %********************
% ZS = linkage(pdist(US,'euclidean'),'ward'); 
% subplot(2,1,1)
% dendrogram(ZS,nbneuronS ); set (gca,'xtick', []);
% title('Dendrogram','fontsize',14); 
% subplot(2,1,2)
% bar(ZS(nbneuronA-1:-1:nbneuronA-20,3));
% title('data number for the last 20 clusters','fontsize',14);
% %print -dpng '/usr/lodyc/tropic/sa/DATA/SMOS/SMOS2016/SVD/figures/KOHONEN_V2018/Kohonen_time_eddy_dendroS_v2018.png' 
% %close
% 
% % on essaie 3 classes
% n=3;
% class_ref_S =cluster(ZS,n);
% som_cplane (sMapS, class_ref_S,1);
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_mapS_v2018.png';
% %close
% % som_bmus va permettre d'associer chaque individu à son référent puis sa classe
% bmusS=som_bmus(sMapS,XS); classS=class_ref_S(bmusS);
% 
% % representation temporelle des classes pour la SSS
% for i=1:length(classS)
%     if classS(i)==1
%            ClassS1(i)=classS(i); ClassS2(i)=NaN; ClassS3(i)=NaN;
%     elseif classS(i)==2
%            ClassS2(i)=classS(i); ClassS3(i)=NaN; ClassS1(i)=NaN;
%     else classS(i)==3
%            ClassS3(i)=classS(i); ClassS2(i)=NaN; ClassS1(i)=NaN;
%     end
% end
% 
% % Trace l occurence des classes de S en fonction du temps
% plot(dates,ClassS1,'c*',dates,ClassS2,'r*',dates,ClassS3,'b*')
% axis([dates(1) dates(nts) 0 n+1]);
% set(gca,'YTick',[0 1 2 3 4]); set(gca,'XTick',[dates(1):365:dates(nts)]); 
% datetick('x','yyyy','keeplimits');
% title('Time serie of SSS classes','fontsize',11);
% xlabel('Time','fontsize',9); ylabel('Class number','fontsize',9);
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_occurenceS_v2018.png';
% %close
% %%
% % reconstruction de la carte S pour la classe 1
% j1=find(ClassS1==1);
% S1=XcarteS(j1,:);
% nbS1=length(j1);
% 
% ipt=0;
% for j=1:nys
% for i=1:nxs
% if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Cs1(i,j,1:nbS1) = S1(1:nbS1, ipt); % si valeur ocean du masque donc 1 on place la valeur de l EOF
% else
% Cs1(i,j,1:nbS1) = NaN;
% end %if
% end %for ix
% end %for iy
% MCs1=nanmean(Cs1,3);
% 
% % reconstruction de la carte S pour la classe 2
% j2=find(ClassS2==2);
% nbS2=length(j2);
% S2=XcarteS(j2,:);
% 
% ipt=0;
% for j=1:nys
% for i=1:nxs
% if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Cs2(i,j,1:nbS2) = S2(1:nbS2, ipt); % si valeur ocean du masque donc 1 on place la valeur de l EOF
% else
% Cs2(i,j,1:nbS2) = NaN;
% end %if
% end %for ix
% end %for iy
% MCs2=nanmean(Cs2,3);
% 
% % reconstruction de la carte S pour la classe 3
% j3=find(ClassS3==3);
% nbS3=length(j3);
% S3=XcarteS(j3,:);
% 
% ipt=0;
% for j=1:nys
% for i=1:nxs
% if (isnan(Masque(i,j))==0 ) ipt = ipt + 1; Cs3(i,j,1:nbS3) = S3(1:nbS3, ipt); % si valeur ocean du masque donc 1 on place la valeur de l EOF
% else
% Cs3(i,j,1:nbS3) = NaN;
% end %if
% end %for ix
% end %for iy
% MCs3=nanmean(Cs3,3);
% %%
% clear j1 j2 j3;
% 
% % affichage des moyennes des 3 classes de la SSS
% set (gca,'fontsize',8)
% subplot(3,2,[1 2])
% contourf(xlos,xlas,MCs1',-1:0.1:0.6)
% axis equal; 
% axis([-70 -30 0 15]);
% colorbar;
% %xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title('SSS mean class 1','fontweight','bold','fontsize',10);
% 
% subplot(3,2,[3 4])
% contourf(xlos,xlas,MCs2',-1:0.1:0.6);
% axis equal; 
% axis([-70 -30 0 15]);
% colorbar;
% %xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title('SSS mean class 2','fontweight','bold','fontsize',10);
% 
% subplot(3,2,[5 6])
% contourf(xlos,xlas,MCs3',-1:0.1:0.6);
% axis equal; 
% axis([-70 -30 0 15]);
% colorbar;
% xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title('SSS mean class 3','fontweight','bold','fontsize',10);
% 
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_meanclassS_v2018.png';
% %close
% %%
% %************************************************************************************************
% %   Les formes fortes
% %************************************************************************************************
% 
% GC=[classT,classA,classS];
% 
% % Identification des formes fortes (cff)
% % Les formes fortes ce sont des 'permutations' des classes, par ex Tclass1 Aclass2 Sclass2 ou Tclass3 Aclass2 Sclass1 etc etc
% cff  = unique(GC,'rows');  % Les classes des formes fortes par ex cff1=1 2 2, cff2=1 2 3  cff3=1 3 1
% nbff = size(cff,1);        % Le nombre de formes fortes
% 
% % Boucle de comptages des formes fortes 
% ind=[]; % sert a reperer les indices des differentes ff pr les entrer ds C
% for i = 1:nbff;  % Pour chaque forme forte :
%     a = ismember(GC, cff(i,:), 'rows'); % On regarde les points des differentes classes qui font partie de la forme forte cad ismember=1
%     b = find(a==1);          % On récupère leurs indices
%     ind=[ind,b'];
%     cardff(i) = length(b);   % Cardinalité de la ième ff 
% end    
% %%
% 
% C1=ind(1:cardff(1)); % A partir des indices cardiff qui donnent le nb de pts de chaque ff on met ds Ci les indices des points de ffi
% compteur=0;
% for i=2:nbff
% compteur=compteur+cardff(i-1);
% eval(['C' num2str(i) ' = ind(compteur+1:compteur+cardff(i))']); % eval permet de faire des matrices de taille variable, ms pas recommande mieux vaut les cellules. Sinon rentrer a la main C1= C2= ....
% 
% %Tab1(i).C=ind(1:cardff(i));
% TaBB(i).C=ind(compteur+1:compteur+cardff(i));
% end
% 
% 
% 
% %%
% 
%  %Tab1=[];
%  compteur=0;
%  clear TEB
%  TEB=[];
%  for i=1:nbff
%      FFI=[];
%      Tabb1(1,i).C=ind(compteur+1:compteur+cardff(i));
% 
%      FFI=[Tabb1(1,i).C' ones(size(Tabb1(1,i).C))'*i];
%      
%      TEB=[TEB; FFI];
%      compteur=compteur+cardff(i);
%  end
%  
%  
%  
%  
% %%
% %***********************    Le profile temporel   ****************************
% Tab=[C1',ones(size(C1))'; % tab concatene le numero de l indice des points temporels avec le numero de la ff
%      C2',ones(size(C2))'*2;
%      C3',ones(size(C3))'*3;
%      C4',ones(size(C4))'*4;
%      C5',ones(size(C5))'*5;
%      C6',ones(size(C6))'*6;
%      C7',ones(size(C7))'*7;
%      C8',ones(size(C8))'*8;
%      C9',ones(size(C9))'*9;
%      C10',ones(size(C10))'*10;
%      C11',ones(size(C11))'*11;
%      C12',ones(size(C12))'*12;
%      C13',ones(size(C13))'*13;
%      C14',ones(size(C14))'*14;
%      C15',ones(size(C15))'*15;
%      C16',ones(size(C16))'*16;
%      C17',ones(size(C17))'*17;
%      C18',ones(size(C18))'*18;
%      C19',ones(size(C19))'*19;
%      C20',ones(size(C20))'*20;
%      C21',ones(size(C21))'*21;
%      ];
%  
% [ind,var]=size(Tab);
%  
% for i=1:ind
% j=find(TEB(:,1)==i);
% Stab(i,:)=TEB(j,:); % chaque Stab(635,2) donne le numero du pas de temps de 1 a 635 associe a celui de la ff
% end
% 
% % On va mettre les donnees dans leur ff respectives S1...S12
% clear S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 
% S1=zeros(2,ind)'; S2=S1; S3=S1; S4=S1; S5=S1; S6=S1; S7=S1; S8=S1; S9 =S1; S10=S1; S11=S1; S12=S1; S13=S1; S14=S1; S15=S1; S16=S1; S17=S1; S18=S1; S19=S1; S20=S1; S21=S1; ;% Initialisation a zero
% 
% 
% %%
% for i=1:nbff
%     SFF(i).C=zeros(2,ind)';
% end
% 
% for i=1:ind 
%         for j=i+1:nbff  
%             if Stab(i,2)==j
%                 SFF(j,j).C(i,:)=Stab(i,:);
%             else
%                 SFF(i,j).C(i,:)=[Stab(i,1),NaN]; 
%         end
%     end
% end
%     
%         
%         
% %%        
% 
% for i=1:size(Stab,1)
%     if Stab(i,2)==1
%            S1(i,:)=Stab(i,:);
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==2
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=Stab(i,:);
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==3
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=Stab(i,:);
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==4
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=Stab(i,:);
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];  
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==5
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=Stab(i,:);
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];  
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==6
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=Stab(i,:);
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];     
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==7
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=Stab(i,:);
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];    
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==8
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=Stab(i,:);
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];    
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==9
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=Stab(i,:);
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];    
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==10
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=Stab(i,:);
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];   
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%      elseif Stab(i,2)==11
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=Stab(i,:);
%            S12(i,:)=[Stab(i,1),NaN];                     
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==12
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=Stab(i,:);  
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==13
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=Stab(i,:);  
%            S14(i,:)=[Stab(i,1),NaN];     
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==14
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN];     
%            S14(i,:)=Stab(i,:);      
%            S15(i,:)=[Stab(i,1),NaN];     
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==15
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=Stab(i,:);      
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==16
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=[Stab(i,1),NaN];
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=Stab(i,:);      ;     
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==17
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=[Stab(i,1),NaN];
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=Stab(i,:);           
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==18
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=[Stab(i,1),NaN];
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN];
%            S18(i,:)=Stab(i,:);     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==19
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=[Stab(i,1),NaN];
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN]; 
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=Stab(i,:);     
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=[Stab(i,1),NaN];     
%     elseif Stab(i,2)==20
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=[Stab(i,1),NaN];
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN]; 
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];     
%            S20(i,:)=Stab(i,:);      
%            S21(i,:)=[Stab(i,1),NaN];     
%     else 
%            S1(i,:)=[Stab(i,1),NaN];
%            S2(i,:)=[Stab(i,1),NaN];
%            S3(i,:)=[Stab(i,1),NaN];
%            S4(i,:)=[Stab(i,1),NaN];
%            S5(i,:)=[Stab(i,1),NaN];
%            S6(i,:)=[Stab(i,1),NaN];
%            S7(i,:)=[Stab(i,1),NaN];
%            S8(i,:)=[Stab(i,1),NaN];
%            S9(i,:)=[Stab(i,1),NaN];
%            S10(i,:)=[Stab(i,1),NaN];
%            S11(i,:)=[Stab(i,1),NaN];
%            S12(i,:)=[Stab(i,1),NaN];
%            S13(i,:)=[Stab(i,1),NaN]; 
%            S14(i,:)=[Stab(i,1),NaN];         
%            S15(i,:)=[Stab(i,1),NaN];
%            S16(i,:)=[Stab(i,1),NaN];     
%            S17(i,:)=[Stab(i,1),NaN]; 
%            S18(i,:)=[Stab(i,1),NaN];     
%            S19(i,:)=[Stab(i,1),NaN];          
%            S20(i,:)=[Stab(i,1),NaN];     
%            S21(i,:)=Stab(i,:);
%     end
% end
% 
% plot(dates,S1(:,2),'c*',dates,S2(:,2),'r*',dates,S3(:,2),'b*',dates,S4(:,2),'y*',dates,S5(:,2),'g*',dates,S6(:,2),'c*',dates,S7(:,2),'m*',dates,S8(:,2),'k*',dates,S9(:,2),'.r',dates,S10(:,2),'.b',dates,S11(:,2),'.g',dates,S12(:,2),'.y',dates,S13(:,2),'c*',dates,S14(:,2),'r*',dates,S15(:,2),'b*',dates,S16(:,2),'y*',dates,S17(:,2),'g*',dates,S18(:,2),'c*',dates,S19(:,2),'m*',dates,S20(:,2),'k*',dates,S21(:,2),'.r')
% axis([dates(1) dates(nts) 0 nbff+1]);
% set(gca,'YTick',[0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21]); set(gca,'XTick',[dates(1):365:dates(nts)]); 
% datetick('x','yyyy','keeplimits'); 
% xlabel('Time','fontsize',9); ylabel('Form number','fontsize',9);
% title('Time evolution of the different eddy strong forms','fontsize',11);
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_strongforms_occurence_v2018.png';
% %%
% % On va chercher quels sont les groupes les plus forts et a quelles dates
% clear l NN dat;
% l=isnan(S1(:,2)); % car la forme est en NaN pour les points ne lui appartenant pas...
% NN=find(l==0);
% dat1=dates(NN);    
% %datestr(dat1)    
% clear l NN 
% 
% l=isnan(S2(:,2));
% NN=find(l==0);
% dat2=dates(NN);    
% %datestr(dat2)    
% clear l NN 
% 
% l=isnan(S3(:,2));
% NN=find(l==0);
% dat3=dates(NN);    
% %datestr(dat3)    
% clear l NN 
% 
% l=isnan(S4(:,2));
% NN=find(l==0);
% dat4=dates(NN);    
% %datestr(dat4)    
% clear l NN 
% 
% l=isnan(S5(:,2));
% NN=find(l==0);
% dat5=dates(NN);    
% %datestr(dat5)    
% clear l NN 
% 
% l=isnan(S6(:,2));
% NN=find(l==0);
% dat6=dates(NN);    
% %datestr(dat6)    
% clear l NN 
% 
% l=isnan(S7(:,2));
% NN=find(l==0);
% dat7=dates(NN);    
% %datestr(dat7)    
% clear l NN
% 
% l=isnan(S8(:,2));
% NN=find(l==0);
% dat8=dates(NN);    
% %datestr(dat8)    
% clear l NN 
% 
% l=isnan(S9(:,2));
% NN=find(l==0);
% dat9=dates(NN);    
% %datestr(dat9)    
% clear l NN 
% 
% l=isnan(S10(:,2));
% NN=find(l==0);
% dat10=dates(NN);    
% %datestr(dat10)    
% clear l NN 
% 
% l=isnan(S11(:,2));
% NN=find(l==0);
% dat11=dates(NN);    
% %datestr(dat11)    
% clear l NN 
% 
% l=isnan(S12(:,2));
% NN=find(l==0);
% dat12=dates(NN);    
% %datestr(dat12)    
% clear l NN 
% 
% l=isnan(S13(:,2));
% NN=find(l==0);
% dat13=dates(NN);    
% %datestr(dat13)    
% clear l NN 
% 
% l=isnan(S14(:,2));
% NN=find(l==0);
% dat14=dates(NN);    
% %datestr(dat14)    
% clear l NN 
% 
% l=isnan(S15(:,2));
% NN=find(l==0);
% dat15=dates(NN);    
% %datestr(dat15)    
% clear l NN 
% 
% l=isnan(S16(:,2));
% NN=find(l==0);
% dat16=dates(NN);    
% %datestr(dat16)    
% clear l NN 
% 
% l=isnan(S17(:,2));
% NN=find(l==0);
% dat17=dates(NN);    
% %datestr(dat17)    
% clear l NN 
% 
% l=isnan(S18(:,2));
% NN=find(l==0);
% dat18=dates(NN);    
% %datestr(dat18)    
% clear l NN 
% 
% l=isnan(S19(:,2));
% NN=find(l==0);
% dat19=dates(NN);    
% %datestr(dat19)    
% clear l NN 
% 
% l=isnan(S20(:,2));
% NN=find(l==0);
% dat20=dates(NN);    
% %datestr(dat20)    
% clear l NN 
% 
% l=isnan(S21(:,2));
% NN=find(l==0);
% dat21=dates(NN);    
% %datestr(dat21)    
% clear l NN 
% 
% dat=[length(dat1) length(dat2) length(dat3) length(dat4) length(dat5) length(dat6) length(dat7) length(dat8) length(dat9) length(dat10) length(dat11) length(dat12) length(dat13) length(dat14) length(dat15) length(dat16) length(dat17) length(dat18) length(dat19) length(dat20) length(dat21)];
% tri=sort(dat,'ascend') %tri par ordre d importance les ff
% for i=1:nbff
% if tri(i)==length(dat1) disp('form number 1')           
% elseif tri(i)==length(dat2) disp('form number 2')
% elseif tri(i)==length(dat3) disp('form number 3')
% elseif tri(i)==length(dat4) disp('form number 4')
% elseif tri(i)==length(dat5) disp('form number 5')
% elseif tri(i)==length(dat6) disp('form number 6')
% elseif tri(i)==length(dat7) disp('form number 7')
% elseif tri(i)==length(dat8) disp('form number 8')
% elseif tri(i)==length(dat9) disp('form number 9')
% elseif tri(i)==length(dat10) disp('form number 10')
% elseif tri(i)==length(dat11) disp('form number 11')
% elseif tri(i)==length(dat12) disp('form number 12')
% elseif tri(i)==length(dat13) disp('form number 13')
% elseif tri(i)==length(dat14) disp('form number 14')
% elseif tri(i)==length(dat15) disp('form number 15')
% elseif tri(i)==length(dat16) disp('form number 16')
% elseif tri(i)==length(dat17) disp('form number 17')
% elseif tri(i)==length(dat18) disp('form number 18')
% elseif tri(i)==length(dat19) disp('form number 19')
% elseif tri(i)==length(dat20) disp('form number 20')
% else disp('form number 21')
% end
% end
% 
% % On entre dans vec la succession des formes pour voir leur evolution temporelle
% vec=S1(:,2);
% u=S2(:,2);
% l=find(isnan(u)==0);vec(l)=2;
% clear u l;
% u=S3(:,2);
% l=find(isnan(u)==0);vec(l)=3;
% clear u l;
% u=S4(:,2);
% l=find(isnan(u)==0);vec(l)=4;
% clear u l;
% u=S5(:,2);
% l=find(isnan(u)==0);vec(l)=5;
% clear u l;
% u=S6(:,2);
% l=find(isnan(u)==0);vec(l)=6;
% clear u l;
% u=S7(:,2);
% l=find(isnan(u)==0);vec(l)=7;
% clear u l;
% u=S8(:,2);
% l=find(isnan(u)==0);vec(l)=8;
% clear u l;
% u=S9(:,2);
% l=find(isnan(u)==0);vec(l)=9;
% clear u l;
% u=S10(:,2);
% l=find(isnan(u)==0);vec(l)=10;
% clear u l;
% u=S11(:,2);
% l=find(isnan(u)==0);vec(l)=11;
% clear u l;
% u=S12(:,2);
% l=find(isnan(u)==0);vec(l)=12;
% clear u l;
% u=S13(:,2);
% l=find(isnan(u)==0);vec(l)=13;
% clear u l;
% u=S14(:,2);
% l=find(isnan(u)==0);vec(l)=14;
% clear u l;
% u=S15(:,2);
% l=find(isnan(u)==0);vec(l)=15;
% clear u l;
% u=S16(:,2);
% l=find(isnan(u)==0);vec(l)=16;
% clear u l;
% u=S17(:,2);
% l=find(isnan(u)==0);vec(l)=17;
% clear u l;
% u=S18(:,2);
% l=find(isnan(u)==0);vec(l)=18;
% clear u l;
% u=S19(:,2);
% l=find(isnan(u)==0);vec(l)=19;
% clear u l;
% u=S20(:,2);
% l=find(isnan(u)==0);vec(l)=20;
% clear u l;
% u=S21(:,2);
% l=find(isnan(u)==0);vec(l)=21;
% clear u l;
% 
% % Trace de l"evolution temporelle des groupes
% plot(dates,vec,dates,S1(:,2),'c*',dates,S2(:,2),'r*',dates,S3(:,2),'b*',dates,S4(:,2),'y*',dates,S5(:,2),'m*',dates,S6(:,2),'g*',dates,S7(:,2),'k*',dates,S8(:,2),'yo',dates,S9(:,2),'bo',dates,S10(:,2),'ro',dates,S11(:,2),'go',dates,S12(:,2),'mo',dates,S13(:,2),'c*',dates,S14(:,2),'r*',dates,S15(:,2),'b*',dates,S16(:,2),'y*',dates,S17(:,2),'m*',dates,S18(:,2),'g*',dates,S19(:,2),'k*',dates,S20(:,2),'yo',dates,S21(:,2))
% axis([dates(1) dates(nts) 0 nbff+1]);
% xlabel('Time','fontsize',9); ylabel('Form number','fontsize',9);
% datetick('x','yyyy','keeplimits');
% title(['Time evolution of the different forms between 2010-2015'],'fontweight','bold','fontsize',11)
% %print -dpng '/net/argos/data/vog/frklod/figures/Kohonen_time_evolution_forms__v2018.png';
% 
% % Trace des formes fortes cff
% % Ici on choisit des formes opposees dans le cycle par ex
% n1=2; form1=eval(['S',num2str(n1)]); 
% n2=5; form2=eval(['S',num2str(n2)]);;
% 
% subplot(14,1,[1 2])
% set(gca,'fontsize',7);
% plot(dates,form1(:,2)/8); % donne la seuqence temporelle
% set(gca,'YTick',[]); datetick('x','yyyy','keeplimits')
% title(['Time evolution strong form ',num2str(n1),'     nb pts = ',num2str(cardff(n1)), '    classes(T,A,S)=(', num2str(cff(n1,1)),',',num2str(cff(n1,2)),',',num2str(cff(n1,3)),')'],'fontweight','bold','fontsize',9)
% 
% subplot(14,1,[4 6])
% set(gca,'fontsize',7);
% atracerT=eval(['MCt',num2str(cff(n1,1))']);
% axis equal; 
% contourf(xlos,xlas,atracerT'); colorbar
% axis([-70 -30 0 15])
% ylabel(' Latitude','fontsize',9); %xlabel('Longitude','fontsize',9
% title(['SST class',num2str(cff(n1,1))],'fontweight','bold','fontsize',9)
% 
% subplot(14,1,[8 10])
% set(gca,'fontsize',7);
% atracerA=eval(['MCa',num2str(cff(n1,2))']);
% axis equal; 
% contourf(xlos,xlas,atracerA'); colorbar
% axis([-70 -30 0 15])
% ylabel(' Latitude','fontsize',9); %xlabel('Longitude','fontsize',9);
% title(['ADT class',num2str(cff(n1,2))],'fontweight','bold','fontsize',9)
% 
% subplot(14,1,[12 14])
% set(gca,'fontsize',7);
% atracerS=eval(['MCs',num2str(cff(n1,3))']);
% axis equal; 
% contourf(xlos,xlas,atracerS',[30:0.5:37]); colorbar
% axis([-70 -30 0 15])
% xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title(['SSS class',num2str(cff(n1,3))],'fontweight','bold','fontsize',9)
% %print ('-dpng',['/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_form',num2str(n1),'_v2018.png']);
% %close
% 
% subplot(14,1,[1 2])
% set(gca,'fontsize',7);
% plot(dates,form2(:,2)/8); % donne la seuqence temporelle
% set(gca,'YTick',[]); datetick('x','yyyy','keeplimits')
% title(['Time evolution strong form ',num2str(n2),'     nb pts = ',num2str(cardff(n2)), '    classes(T,A,S)=(', num2str(cff(n2,1)),',',num2str(cff(n2,2)),',',num2str(cff(n2,3)),')'],'fontweight','bold','fontsize',9)
% 
% subplot(14,1,[4 6])
% set(gca,'fontsize',7);
% atracerT=eval(['MCt',num2str(cff(n2,1))']);
% axis equal; 
% contourf(xlos,xlas,atracerT'); colorbar
% axis([-70 -30 0 15])
% ylabel(' Latitude','fontsize',9); %xlabel('Longitude','fontsize',9
% title(['SST class',num2str(cff(n2,1))],'fontweight','bold','fontsize',9)
% 
% subplot(14,1,[8 10])
% set(gca,'fontsize',7);
% atracerA=eval(['MCa',num2str(cff(n2,2))']);
% axis equal; 
% contourf(xlos,xlas,atracerA'); colorbar
% axis([-70 -30 0 15])
% ylabel(' Latitude','fontsize',9); %xlabel('Longitude','fontsize',9);
% title(['ADT class',num2str(cff(n2,2))],'fontweight','bold','fontsize',9)
% 
% subplot(14,1,[12 14])
% set(gca,'fontsize',7);
% atracerS=eval(['MCs',num2str(cff(n2,3))']);
% axis equal; 
% contourf(xlos,xlas,atracerS',[30:0.5:37]); colorbar
% axis([-70 -30 0 15])
% xlabel('Longitude','fontsize',9);ylabel(' Latitude','fontsize',9); 
% title(['SSS class',num2str(cff(n2,3))],'fontweight','bold','fontsize',9)
% %print ('-dpng',['/net/argos/data/vog/frklod/figures/Kohonen_time_eddy_form',num2str(n2),'_v2018.png']);
% %close
% 
