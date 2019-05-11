%Pour Statitiques descriptives Finales
clear all ; close all; clc ;

load('StatDescriptivesClasses.mat')
nt = 635; 
XSS1=reshape(class1_sss_all,[nxs*nys,nt]);



XS=XSS1(ij_masquet_sst_bon,jour); 
XS=XS';


%%% Base de test 
load('dates.mat', 'dates')

tout = [1:635];

ValidationTest = setdiff(tout,jourTests);

% Test des jours pour coonfirmation
jourValidationTest =datestr(dates(ValidationTest));


load XSS.mat
load XTT.mat
load XAA.mat
load ij_masquet_sst_bon.mat

 Tot_SSS=[]; Tot_SST=[];Tot_ADT=[];
for jour = 1:length(ValidationTest)
    XT=XTT(ij_masquet_sst_bon,ValidationTest(jour));
    XT=XT';
    
    XS=XSS(ij_masquet_sst_bon,ValidationTest(jour));
    XS=XS';
    
    XA=XAA(ij_masquet_sst_bon,ValidationTest(jour));
    XA=XA';
       
    for iDS = 1
        Tot_SSS=[Tot_SSS;XS(iDS,:)'];
    end
    
    for iDS = 1
        Tot_SST=[Tot_SST;XT(iDS,:)'];
    end
    
    
    for iDS = 1
        Tot_ADT=[Tot_ADT;XA(iDS,:)'];
    end  
    
end

Xvt = [Tot_SSS Tot_SST Tot_ADT];

save baseValidationTest Xvt ValidationTest jourValidationTest



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot(class1_sss,'DisplayName','class1_sss');hold on;
% plot(class2_sss,'DisplayName','class2_sss');
% plot(class3_sss,'DisplayName','class3_sss');
% hold off;
Tabclass1 = [class1_sss class1_sst class1_adt];

TabSSS = [class1_sss class2_sss class3_sss];
StatSSS = [round(nanmean(TabSSS)',2) round(nanmin(TabSSS)',2) round(nanmax(TabSSS)',2) round(nanstd(TabSSS)',2)];

TabSST = [class1_sst class2_sst class3_sst];
StatSST = [round(nanmean(TabSST)',2) round(nanmin(TabSST)',2) round(nanmax(TabSST)',2) round(nanstd(TabSST)',2)];

TabADT = [class1_adt class2_adt class3_adt];
StatADT = [round(nanmean(TabADT)',2) round(nanmin(TabADT)',2) round(nanmax(TabADT)',2) round(nanstd(TabADT)',2)];
