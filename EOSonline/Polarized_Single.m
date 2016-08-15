%%
addpath('../Library');
warning('off','all');
%%
filename='/Users/Zhenjie/Data/2016-07-25/07-25-2016_19_27_25_TopB.fits';
[Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EF]=EOS_Online( filename ,...
    'ShowOutline',1,'KappaMode',5,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',2.1,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',295,'IfBin',0,'BinGridSize',160,'pixellength',0.7e-6,'SM',30);

%%
d=figure();
scatter(Ptsel,Ktsel);
hold on
xlim([0,6])
[ KappaTildeP, PTildeP, ~, ~ ] = IdealFermiEOS( 1.1, 2, 100 );
ylim([0,4])
plot(PTildeP,KappaTildeP);
hold off

%'TailRange',[65,360],