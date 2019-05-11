function CcDon1 = denorm_min_max(Don_cr, mini, maxi);



%% Critere de d??normalisation de "centre-reduction": X=Xn.*e+m
% m: la moyenne de X
% e: Ecart type de X

%vec_cr= CcDon_cr;
eval(['vec_cr=Don_cr;']);
%longueur = length(vect);
eval(['minim=mini;']);
eval(['maxim=maxi;']);
vecta = vec_cr*(maxim-minim)+ minim;
CcDon1 = vecta;
