%%
addpath('../Library');
warning('off','all');
%%
filename='J:\Elder Backup Raw Images\2016\2016-08\2016-08-15\08-15-2016_17_36_53_TopB.fits';
[Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFsort,P,zcor,Vsel]=EOS_Online( filename,...
    'ShowOutline',1,'KappaMode',5,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',1.2,'kmin',0.,...
    'Fudge',1.7,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',inf,'IfBin',1,'BinGridSize',160,'SM',4,'CleanImage',1,'pixellength',0.7e-6);
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