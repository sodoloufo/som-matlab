function CcDonn1 = normal_min_max(CcDon_test, mini, maxi);




%% Critere de centre reduction: Xn=(X-m)/e
% m: la moyenne de X
% e: Ecart type de X
%% Cette fonction produit le yyy_don_cr en fonction du vecteur yyy_don


eval(['vect=CcDon_test;']);
eval(['minim=mini;']);
eval(['maxim=maxi;']);
vec_cr=(vect - minim)/(maxim-minim);
CcDonn1=vec_cr;





