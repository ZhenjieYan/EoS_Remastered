function [muTilde_ouput,TTilde_ouput,impurity_ouput] = ...
    MinorityChemicalPotentialVirial( Select_TTilde1, Select_impurity1 )

% Select_TTilde1:T/TF1
% Select_impurity1: n1/n2
% muTilde_ouput: mu1/TF2;
% TTilde_ouput: T/TF1
% impurity_ouput: 

DictionaryFolder='/Users/Zhenjie/Data/SpinImbalanceVirial/';
load([DictionaryFolder,'SpinImbalanceVirial.mat']);

DTTilde1 = abs(TTilde1 - Select_TTilde1);

Dimpurity1 = abs(impurity1 - Select_impurity1);

Dboth = DTTilde1 + Dimpurity1;

[MinBoth,MinVecIndex]= min(Dboth(:));
[Min_row,Min_col] = ind2sub(size(Dboth),MinVecIndex);

muTilde_ouput = muTilde1(Min_row,Min_col);
TTilde_ouput = TTilde1(Min_row,Min_col);
impurity_ouput = impurity1(Min_row,Min_col);

end
