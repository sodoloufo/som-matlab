class_dat=class_data(1:10:end,:); % classe de donnees
bmus_1=bmus_smap3(find(class_dat==1)); % Pour trouver tous les bmus correspondant à la classe 1
val_b1=[];
for i=1:length(bmus_1)                 % Pour trouver les neurones de la carte smap3 correspondant à bmus de la classe 1
    val_b=smap3.codebook(bmus_1(i),:)
    val_b1=[val_b1;val_b];
end
mean_b1=mean(val_b1);                  % moyenne de la sst et de la chloro