%% Preloading
%Define the physical constant
mli=9.988346*10^-27;  %kg
hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant
omega=23.9*2*pi; %in rad/s
pixellength=0.7*10^-6; %in m
sigma0=0.215*10^-12/2; %in m^2
%load all the functions
addpath('../Library');
filefolder='/Users/Zhenjie/Data/2016-07-25/';
fileS1='07-25-2016_19_27_25_TopA.fits';
fileS2='07-25-2016_19_27_25_TopB.fits';
Fudge=2.2;
%% Get the profile for all of them

[Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( [filefolder,fileS1],'ROI1',[100,200,700,1375],...
    'ROI2',[300,642,712,833],'TailRange',[350,1100],'ShowOutline',0,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',295,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength);
nS1=nsort/1e18;
VS1=Vsort/hh;
Z0S1=zcor.z0*pixellength/1e-6;
ZS1=Zsort/1e-6;

%% 'ROI2',[300,515,712,930],

[Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( [filefolder,fileS2],'ROI1',[100,200,700,1375],...
    'ROI2',[300,617,712,837],'TailRange',[350,1100],'ShowOutline',0,'KappaMode',5,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',295,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',50);
nS2=nsort/1e18;
VS2=Vsort/hh;
Z0S2=zcor.z0*pixellength/1e-6;
ZS2=Zsort/1e-6;

%%
plot(ZS1,nS1,'b.',ZS2,nS2,'r.')
%%
plot(VS1,nS1,'b.',VS2,nS2,'r.');xlim([0,10000])
%%
Vgrid=VS1;
nS1grid=nS1;
nS2grid=interp1(VS2,nS2,Vgrid);

plot(Vgrid,nS1grid,'b.',Vgrid,nS2grid,'r.');xlim([0,10000])
%%
Vth=2000;
mask=Vgrid>Vth;

Vfit=Vgrid(mask);
nS1fit=nS1grid(mask);

P=IdealGasFit( Vfit,nS1fit,[5000,500,1] );
%%
nS1Ideal=IdealGasnV( P(1:2),Vgrid );
plot(Vgrid,nS1grid,'b.',Vgrid,nS1Ideal+0.35*nS2grid,'r.');xlim([0,10000])

