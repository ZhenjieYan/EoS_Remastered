%%
addpath('../Library');
warning('off','all');
%%
filename='/Volumes/Raw Data/Images/2016/2016-10/2016-10-20/10-20-2016_23_40_30_TopA.fits';
[Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EF]=EOS_Online( filename ,...
    'ShowOutline',1,'KappaMode',5,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',1.55,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',inf,'IfBin',1,'BinGridSize',160,'pixellength',0.7e-6,'SM',4,'CleanImage',1,'BoxShape','Square');

