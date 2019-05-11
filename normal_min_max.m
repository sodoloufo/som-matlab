function [CcDonn1, mini, maxi] = normal_min_max(CcDon);




%% Critere de centre reduction: Xn=(X-m)/e
% m: la moyenne de X
% e: Ecart type de X
%% Cette fonction produit le yyy_don_cr en fonction du vecteur yyy_don


eval(['vect=CcDon;']);
eval(['mini=nanmin(reshape(CcDon,[],1));']);
eval(['maxi=nanmax(reshape(CcDon,[],1));']);
vec_cr=(vect - mini)/(maxi-mini);
CcDonn1=vec_cr;




