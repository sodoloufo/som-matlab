clear all ; close all; clc ;
% addpath ('som');
% addpath ('som\dijkstradir');
addpath('base_traite');

%load smap_class.mat
load Masque.mat
load nys.mat
load nxs.mat
load xlos.mat
load xlas.mat
% load mini.mat;
% load maxi.mat;
load dates.mat
% load results_app2_40_30.mat
load XSS.mat
load XTT.mat
load XAA.mat
load ij_masquet_sst_bon.mat

% Wei = sMap.codebook;

% formatIn = 'dd-mm-yyyy';
Dim=[154 76];
CLimSSS=[25 40];
CLimSST=[25 30];
CLimADT=[0.2 0.8];

% % Traitement pour climatologie

load('DatesEtRepartMoisSaisons.mat')
% % % %%%%% -----------------Bases mensuelles SST ------------------------

nomJour = 'Janvier';
XTM1=XTT(ij_masquet_sst_bon,indexM1);
MoyXTM1=mean(XTM1,2);
% MoyXTM1 = MoyXTM1';
plot7_imageMois('Mesure moyenne SST',MoyXTM1',ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Fevrier';
XTM2=XTT(ij_masquet_sst_bon,indexM2);
MoyXTM2=mean(XTM2,2);
MoyXTM2 = MoyXTM2';
plot7_imageMois('Mesure moyenne SST',MoyXTM2(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Mars';
XTM3=XTT(ij_masquet_sst_bon,indexM3);
MoyXTM3=mean(XTM3,2);
MoyXTM3 = MoyXTM3';
plot7_imageMois('Mesure moyenne SST',MoyXTM3(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Avril';
XTM4=XTT(ij_masquet_sst_bon,indexM4);
MoyXTM4=mean(XTM4,2);
MoyXTM4 = MoyXTM4';
plot7_imageMois('Mesure moyenne SST',MoyXTM4(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Mai';
XTM5=XTT(ij_masquet_sst_bon,indexM5);
MoyXTM5=mean(XTM5,2);
MoyXTM5 = MoyXTM5';
plot7_imageMois('Mesure moyenne SST',MoyXTM5(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Juin';
XTM6=XTT(ij_masquet_sst_bon,indexM6);
MoyXTM6=mean(XTM6,2);
MoyXTM6 = MoyXTM6';
plot7_imageMois('Mesure moyenne SST',MoyXTM6(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);


nomJour = 'Juillet';
XTM7=XTT(ij_masquet_sst_bon,indexM7);
MoyXTM7=mean(XTM7,2);
MoyXTM7 = MoyXTM7';
plot7_imageMois('Mesure moyenne SST',MoyXTM7(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Aout';
XTM8=XTT(ij_masquet_sst_bon,indexM8);
MoyXTM8=mean(XTM8,2);
MoyXTM8 = MoyXTM8';
plot7_imageMois('Mesure moyenne SST',MoyXTM8(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Septembre';
XTM9=XTT(ij_masquet_sst_bon,indexM9);
MoyXTM9=mean(XTM9,2);
MoyXTM9 = MoyXTM9';
plot7_imageMois('Mesure moyenne SST',MoyXTM9(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Octobre';
XTM10=XTT(ij_masquet_sst_bon,indexM10);
MoyXTM10=mean(XTM10,2);
MoyXTM10 = MoyXTM10';
plot7_imageMois('Mesure moyenne SST',MoyXTM10(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Novembre';
XTM11=XTT(ij_masquet_sst_bon,indexM11);
MoyXTM11=mean(XTM11,2);
MoyXTM11 = MoyXTM11';
plot7_imageMois('Mesure moyenne SST',MoyXTM11(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Decembre';
XTM12=XTT(ij_masquet_sst_bon,indexM12);
MoyXTM12=mean(XTM12,2);
MoyXTM12 = MoyXTM12';
plot7_imageMois('Mesure moyenne SST',MoyXTM12(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

% % Par saison

nomJour = 'Hiver';
XTHiver=XTT(ij_masquet_sst_bon,indexHiver);
MoyXTHiver=mean(XTHiver,2);
MoyXTHiver = MoyXTHiver';
plot7_imageMois('Mesure moyenne SST',MoyXTHiver(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Ete';
XTEte=XTT(ij_masquet_sst_bon,indexEte);
MoyXTEte=mean(XTEte,2);
MoyXTEte = MoyXTEte';
plot7_imageMois('Mesure moyenne SST',MoyXTEte(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Automne';
XTAutom=XTT(ij_masquet_sst_bon,indexAutom);
MoyXTAutom=mean(XTAutom,2);
MoyXTAutom = MoyXTAutom';
plot7_imageMois('Mesure moyenne SST',MoyXTAutom(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);

nomJour = 'Printemps';
XTPrint=XTT(ij_masquet_sst_bon,indexPrint);
MoyXTPrint=mean(XTPrint,2);
MoyXTPrint = MoyXTPrint';
plot7_imageMois('Mesure moyenne SST',MoyXTPrint(1,:),ij_masquet_sst_bon,CLimSST,nomJour,xlos,xlas,Dim);



% % % %%%%% -----------------Bases mensuelles SSS ------------------------

nomJour = 'Janvier';
XSM1=XSS(ij_masquet_sst_bon,indexM1);
MoyXSM1=mean(XSM1,2);
% MoyXSM1 = MoyXSM1';
plot7_imageMois('Mesure moyenne SSS',MoyXSM1',ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Fevrier';
XSM2=XSS(ij_masquet_sst_bon,indexM2);
MoyXSM2=mean(XSM2,2);
MoyXSM2 = MoyXSM2';
plot7_imageMois('Mesure moyenne SSS',MoyXSM2(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Mars';
XSM3=XSS(ij_masquet_sst_bon,indexM3);
MoyXSM3=mean(XSM3,2);
MoyXSM3 = MoyXSM3';
plot7_imageMois('Mesure moyenne SSS',MoyXSM3(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Avril';
XSM4=XSS(ij_masquet_sst_bon,indexM4);
MoyXSM4=mean(XSM4,2);
MoyXSM4 = MoyXSM4';
plot7_imageMois('Mesure moyenne SSS',MoyXSM4(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Mai';
XSM5=XSS(ij_masquet_sst_bon,indexM5);
MoyXSM5=mean(XSM5,2);
MoyXSM5 = MoyXSM5';
plot7_imageMois('Mesure moyenne SSS',MoyXSM5(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Juin';
XSM6=XSS(ij_masquet_sst_bon,indexM6);
MoyXSM6=mean(XSM6,2);
MoyXSM6 = MoyXSM6';
plot7_imageMois('Mesure moyenne SSS',MoyXSM6(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);


nomJour = 'Juillet';
XSM7=XSS(ij_masquet_sst_bon,indexM7);
MoyXSM7=mean(XSM7,2);
MoyXSM7 = MoyXSM7';
plot7_imageMois('Mesure moyenne SSS',MoyXSM7(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Aout';
XSM8=XSS(ij_masquet_sst_bon,indexM8);
MoyXSM8=mean(XSM8,2);
MoyXSM8 = MoyXSM8';
plot7_imageMois('Mesure moyenne SSS',MoyXSM8(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Septembre';
XSM9=XSS(ij_masquet_sst_bon,indexM9);
MoyXSM9=mean(XSM9,2);
MoyXSM9 = MoyXSM9';
plot7_imageMois('Mesure moyenne SSS',MoyXSM9(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Octobre';
XSM10=XSS(ij_masquet_sst_bon,indexM10);
MoyXSM10=mean(XSM10,2);
MoyXSM10 = MoyXSM10';
plot7_imageMois('Mesure moyenne SSS',MoyXSM10(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Novembre';
XSM11=XSS(ij_masquet_sst_bon,indexM11);
MoyXSM11=mean(XSM11,2);
MoyXSM11 = MoyXSM11';
plot7_imageMois('Mesure moyenne SSS',MoyXSM11(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Decembre';
XSM12=XSS(ij_masquet_sst_bon,indexM12);
MoyXSM12=mean(XSM12,2);
MoyXSM12 = MoyXSM12';
plot7_imageMois('Mesure moyenne SSS',MoyXSM12(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

% % Par saison

nomJour = 'Hiver';
XSHiver=XSS(ij_masquet_sst_bon,indexHiver);
MoyXSHiver=mean(XSHiver,2);
MoyXSHiver = MoyXSHiver';
plot7_imageMois('Mesure moyenne SSS',MoyXSHiver(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Ete';
XSEte=XSS(ij_masquet_sst_bon,indexEte);
MoyXSEte=mean(XSEte,2);
MoyXSEte = MoyXSEte';
plot7_imageMois('Mesure moyenne SSS',MoyXSEte(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Automne';
XSAutom=XSS(ij_masquet_sst_bon,indexAutom);
MoyXSAutom=mean(XSAutom,2);
MoyXSAutom = MoyXSAutom';
plot7_imageMois('Mesure moyenne SSS',MoyXSAutom(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

nomJour = 'Printemps';
XSPrint=XSS(ij_masquet_sst_bon,indexPrint);
MoyXSPrint=mean(XSPrint,2);
MoyXSPrint = MoyXSPrint';
plot7_imageMois('Mesure moyenne SSS',MoyXSPrint(1,:),ij_masquet_sst_bon,CLimSSS,nomJour,xlos,xlas,Dim);

% % % %%%%% -----------------Bases mensuelles ADT ------------------------

nomJour = 'Janvier';
XAM1=XAA(ij_masquet_sst_bon,indexM1);
MoyXAM1=mean(XAM1,2);
MoyXAM1 = MoyXAM1';
plot7_imageMois('Mesure moyenne ADT',MoyXAM1(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Fevrier';
XAM2=XAA(ij_masquet_sst_bon,indexM2);
MoyXAM2=mean(XAM2,2);
MoyXAM2 = MoyXAM2';
plot7_imageMois('Mesure moyenne ADT',MoyXAM2(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Mars';
XAM3=XAA(ij_masquet_sst_bon,indexM3);
MoyXAM3=mean(XAM3,2);
MoyXAM3 = MoyXAM3';
plot7_imageMois('Mesure moyenne ADT',MoyXAM3(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Avril';
XAM4=XAA(ij_masquet_sst_bon,indexM4);
MoyXAM4=mean(XAM4,2);
MoyXAM4 = MoyXAM4';
plot7_imageMois('Mesure moyenne ADT',MoyXAM4(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Mai';
XAM5=XAA(ij_masquet_sst_bon,indexM5);
MoyXAM5=mean(XAM5,2);
MoyXAM5 = MoyXAM5';
plot7_imageMois('Mesure moyenne ADT',MoyXAM5(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Juin';
XAM6=XAA(ij_masquet_sst_bon,indexM6);
MoyXAM6=mean(XAM6,2);
MoyXAM6 = MoyXAM6';
plot7_imageMois('Mesure moyenne ADT',MoyXAM6(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);


nomJour = 'Juillet';
XAM7=XAA(ij_masquet_sst_bon,indexM7);
MoyXAM7=mean(XAM7,2);
MoyXAM7 = MoyXAM7';
plot7_imageMois('Mesure moyenne ADT',MoyXAM7(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Aout';
XAM8=XAA(ij_masquet_sst_bon,indexM8);
MoyXAM8=mean(XAM8,2);
MoyXAM8 = MoyXAM8';
plot7_imageMois('Mesure moyenne ADT',MoyXAM8(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Septembre';
XAM9=XAA(ij_masquet_sst_bon,indexM9);
MoyXAM9=mean(XAM9,2);
MoyXAM9 = MoyXAM9';
plot7_imageMois('Mesure moyenne ADT',MoyXAM9(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Octobre';
XAM10=XAA(ij_masquet_sst_bon,indexM10);
MoyXAM10=mean(XAM10,2);
MoyXAM10 = MoyXAM10';
plot7_imageMois('Mesure moyenne ADT',MoyXAM10(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Novembre';
XAM11=XAA(ij_masquet_sst_bon,indexM11);
MoyXAM11=mean(XAM11,2);
MoyXAM11 = MoyXAM11';
plot7_imageMois('Mesure moyenne ADT',MoyXAM11(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Decembre';
XAM12=XAA(ij_masquet_sst_bon,indexM12);
MoyXAM12=mean(XAM12,2);
MoyXAM12 = MoyXAM12';
plot7_imageMois('Mesure moyenne ADT',MoyXAM12(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

% % Par saison

nomJour = 'Hiver';
XAHiver=XAA(ij_masquet_sst_bon,indexHiver);
MoyXAHiver=mean(XAHiver,2);
MoyXAHiver = MoyXAHiver';
plot7_imageMois('Mesure moyenne ADT',MoyXAHiver(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Ete';
XAEte=XAA(ij_masquet_sst_bon,indexEte);
MoyXAEte=mean(XAEte,2);
MoyXAEte = MoyXAEte';
plot7_imageMois('Mesure moyenne ADT',MoyXAEte(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Automne';
XAAutom=XAA(ij_masquet_sst_bon,indexAutom);
MoyXAAutom=mean(XAAutom,2);
MoyXAAutom = MoyXAAutom';
plot7_imageMois('Mesure moyenne ADT',MoyXAAutom(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);

nomJour = 'Printemps';
XAPrint=XAA(ij_masquet_sst_bon,indexPrint);
MoyXAPrint=mean(XAPrint,2);
MoyXAPrint = MoyXAPrint';
plot7_imageMois('Mesure moyenne ADT',MoyXAPrint(1,:),ij_masquet_sst_bon,CLimADT,nomJour,xlos,xlas,Dim);


% % Moyennes mensuelles


MoyADT=[];
MoyADT=[mean(MoyXAM12);mean(MoyXAM1);mean(MoyXAM2);mean(MoyXAM3);...
    mean(MoyXAM4);mean(MoyXAM5);mean(MoyXAM6);mean(MoyXAM7);...
    mean(MoyXAM8);mean(MoyXAM9);mean(MoyXAM10);mean(MoyXAM11)];

MoySST=[];
MoySST=[mean(MoyXTM12);mean(MoyXTM1);mean(MoyXTM2);mean(MoyXTM3);...
    mean(MoyXTM4);mean(MoyXTM5);mean(MoyXTM6);mean(MoyXTM7);...
    mean(MoyXTM8);mean(MoyXTM9);mean(MoyXTM10);mean(MoyXTM11)];

MoySSS=[];
MoySSS=[mean(MoyXSM12);mean(MoyXSM1);mean(MoyXSM2);mean(MoyXSM3);...
    mean(MoyXSM4);mean(MoyXSM5);mean(MoyXSM6);mean(MoyXSM7);...
    mean(MoyXSM8);mean(MoyXSM9);mean(MoyXSM10);mean(MoyXSM11)];

%%%%%%%
StatMoy = [round(MoySSS,2) round(MoySST,2) round(MoyADT,2)];

StatSaisonsXA = [round(mean(MoyXAHiver),2) round(mean(MoyXAPrint),2) round(mean(MoyXAEte),2) round(mean(MoyXAAutom),2)];

StatSaisonsXS = [round(mean(MoyXSHiver),2) round(mean(MoyXSPrint),2) round(mean(MoyXSEte),2) round(mean(MoyXSAutom),2)];
StatSaisonsXT = [round(mean(MoyXTHiver),2) round(mean(MoyXTPrint),2) round(mean(MoyXTEte),2) round(mean(MoyXTAutom),2)];

y1 = MoySSS';
y2 = MoySST';
y3 = MoyADT';

figure
plot(x,y1,'g',x,y2,'b--o',x,y3,'c*')

figure
%plot(x,y1,'g')
plot(x,y1,'-o','MarkerIndices',1:1:length(y1))

figure
plot(x,y1,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])



%%%%
x = 1:1:12;
figure
ax1 = subplot(3,1,1); % top subplot
plot(ax1,x,y1,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title(ax1,'Moyennes mensuelles SSS','fontsize',18)
ylabel(ax1,'Valeurs (psu)')

ax2 = subplot(3,1,2); % bottom subplot
plot(ax2,x,y2,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title(ax2,'Moyennes mensuelles SST','fontsize',18)
ylabel(ax2,'Valeurs (Degre)')

ax3 = subplot(3,1,3); % bottom subplot
plot(ax3,x,y3,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5])
title(ax3,'Moyennes mensuelles ADT','fontsize',18)
ylabel(ax3,'Valeurs (m)')
xlabel('Mois (de Decembre a Novembre)','FontSize',18,'FontWeight','bold')

FileName=[pwd, '/figuresmensuelles/statmensuelles'];
fig = gcf; fig.PaperUnits = 'inches';
fig.PaperPosition = [0 2.5 8.5 3.2];
print(FileName,'-dpng','-r0')
% close

load('BasesMoyennesTempo.mat')
%%coeff = pca(StatMoy);