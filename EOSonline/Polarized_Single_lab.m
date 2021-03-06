%%
addpath('../Library');
warning('off','all');
%%
filename='J:\Elder Backup Raw Images\2016\2016-09\2016-09-01\09-01-2016_19_20_19_TopB.fits';
rawimg=fitsread(filename);
rawBin=ImgBin(rawimg,3);
[Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFsort,P,zcor,Vsel]=EOS_Online( rawBin,...
    'ShowOutline',1,'KappaMode',5,'PolyOrder',10,'VrangeFactor',5,'IfHalf',1,'kmax',1.2,'kmin',0.,...
    'Fudge',3.5,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',inf,'IfBin',1,'BinGridSize',200,'SM',8,'CleanImage',1,'pixellength',0.7e-6,'IfLookUpTable',1);
nS1=nsort;
VS1=Vsort;
ZS1=Zsort;
%%
d=figure();
scatter(Ptsel,Ktsel);
hold on
xlim([0,6])
[ KappaTildeP, PTildeP, ~, ~ ] = IdealFermiEOS();
ylim([0,4])
plot(PTildeP,KappaTildeP);
hold off

%'TailRange',[65,360],