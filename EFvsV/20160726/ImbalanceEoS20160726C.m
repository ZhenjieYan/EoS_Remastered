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
Fudge=1.4344;
kB=1.380e-23;
Nsat=120;
load('Data/2016-07-26C.mat')
warning ('off','all')
%% Get the profile for all of them
nS1list={};
ZS1list={};
Z0S1list=[];
for i=1:length(imglistS1Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( imglistS1Final{i},'ROI1',[100,300,700,1175],...
    'ROI2',[300,642,712,833],'TailRange',[350,1100],'ShowOutline',0,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,'IfCleanImage',1);
    nS1list=[nS1list;nsort/1e18];
    Z0S1=zcor.z0*pixellength/1e-6;
    Z0S1list=[Z0S1list;Z0S1];
    ZS1=Zsort/1e-6+Z0S1;
    ZS1list=[ZS1list;ZS1];
end
%% 'ROI2',[300,515,712,930],
nS2list={};
ZS2list={};
Z0S2list=[];
for i=1:length(imglistS2Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( imglistS2Final{i},'ROI1',[100,300,700,1175],...
    'ROI2',[300,642,712,833],'TailRange',[350,1100],'ShowOutline',0,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,'IfCleanImage',1);
    nS2list=[nS2list;nsort/1e18];
    Z0S2=zcor.z0*pixellength/1e-6;
    Z0S2list=[Z0S2list;Z0S2];
    ZS2=Zsort/1e-6+Z0S2;
    ZS2list=[ZS2list;ZS2];
end
%%
nS1=[];
ZS1=[];
for i=1:length(imglistS1Final)
    nS1=[nS1;nS1list{i}];
    ZS1=[ZS1;ZS1list{i}-mean(Z0S1list)];
end
nS2=[];
ZS2=[];
for i=1:length(imglistS2Final)
    nS2=[nS2;nS2list{i}];
    ZS2=[ZS2;ZS2list{i}-mean(Z0S2list)];
end

UnbinnedProfile.nS1=nS1;
UnbinnedProfile.nS2=nS2;
UnbinnedProfile.ZS1=ZS1;
UnbinnedProfile.ZS2=ZS2;

%% Bin the profile with Z
Nbin=250;
ZGrid=linspace(-300,300,Nbin+1);
[ZS1BinZ,nS1Bin]=BinGrid(ZS1,nS1,ZGrid,0);
[ZS2BinZ,nS2Bin]=BinGrid(ZS2,nS2,ZGrid,0);
%nS1Bin=TailTailor(nS1Bin,ZS1BinZ,-170,170);

plot(ZS1BinZ,nS1Bin,'b',ZS2BinZ,nS2Bin,'r')

ZBinZ=ZS1BinZ;
nS1BinZ=nS1Bin;
nS2BinZ=interp1(ZS2BinZ,nS2Bin,ZBinZ,'spline');
VBinZ=0.5*mli*omega^2*(ZBinZ*1e-6).^2/hh;
plot(ZBinZ,nS1BinZ,'b.',ZBinZ,nS2BinZ,'r.');
plot(VBinZ,nS1BinZ,'b.',VBinZ,nS2BinZ,'r.');
ntotBinZ=nS1BinZ+nS2BinZ;
EFS1BinZ=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS1BinZ)*1e18).^(2/3).*sign(nS1BinZ)/hh;
EFS2BinZ=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS2BinZ)*1e18).^(2/3).*sign(nS1BinZ)/hh;

Zprofile.Z=ZBinZ;
Zprofile.nS1=nS1BinZ;
Zprofile.nS2=nS2BinZ;
Zprofile.V=VBinZ;
Zprofile.EFS1=EFS1BinZ;
Zprofile.EFS2=EFS2BinZ;
%% Bin the profile with V
VS1=0.5*mli*omega^2*(ZS1*1e-6).^2/hh;
VS2=0.5*mli*omega^2*(ZS2*1e-6).^2/hh;
EFS1=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS1)*1e18).^(2/3).*sign(nS1)/hh;
EFS2=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS2)*1e18).^(2/3).*sign(nS1)/hh;

Nbin=200;
VGrid=linspace(0,1.2e4,Nbin+1);
[VS1BinV,EFS1BinV]=BinGrid(VS1,EFS1,VGrid,0);
mask=isnan(VS1BinV);
VS1BinV(mask)=[];
[VS2BinV,EFS2BinV1]=BinGrid(VS2,EFS2,VGrid,0);
mask=isnan(VS2BinV);
VS2BinV(mask)=[];

VBinV=VS1BinV;
%nS1BG=mean(nS1BinV(VBinV>8000));
EFS2BinV=interp1(VS2BinV,EFS2BinV1,VBinV);
nS1BinV=(2*mli*abs(EFS1BinV)*hh/hbar^2).^(3/2)/(6*pi^2)/1e18.*sign(EFS1BinV);
nS2BinV=(2*mli*abs(EFS2BinV)*hh/hbar^2).^(3/2)/(6*pi^2)/1e18.*sign(EFS2BinV);

%nS2BinV=nS2BinV-nS2BG;


ZBinV=sqrt(2*VBinV*hh/(mli*omega^2))/1e-6;


plot(VBinV,EFS1BinV,'r.',VBinV,EFS2BinV,'b.');
%plot(VBinV,nS1BinV,'r.',VBinV,nS2BinV,'b.');

Vprofile.V=VBinV;
Vprofile.nS1=nS1BinV;
Vprofile.nS2=nS2BinV;
Vprofile.EFS1=EFS1BinV;
Vprofile.EFS2=EFS2BinV;
Vprofile.Z=ZBinV;
%%
Vth=1500;
mask=VBinZ>Vth;

Vfit=VBinZ(mask);
nS1fit=nS1BinZ(mask);

P=IdealGasFit( Vfit,nS1fit,[5000,500,1] );
%%
mu_trap_S1=P(1)
T_trap=P(2)
%%
% Get the EoS of ideal fermi gas
[ KappaTildeT, PTildeT, TTildeT, CVTildeT , beta_mu_T ,ZTildeT] = IdealFermiEOS( );
mu_EF_T=beta_mu_T.*TTildeT;
mu_localZ=mu_trap_S1-VBinZ;
% get the kappa from the known T and mu
beta_mu_localZ=mu_localZ/T_trap;
KappaFitZ=interp1(beta_mu_T,KappaTildeT,beta_mu_localZ,'spline');
% get the kappa from image profile
[KappaTildeS1BinV,~]=FiniteD( VBinV,VBinV*0,EFS1BinV,EFS1BinV*0,10);
KappaTildeS1BinV=-KappaTildeS1BinV;
plot(VBinV,KappaTildeS1BinV,'bo');
hold on
plot(VBinZ,KappaFitZ,'k-')
hold off
%%

[KappaTildeS2BinV,~]=FiniteD( VBinV,VBinV*0,EFS2BinV,EFS2BinV*0,2);
KappaTildeS2BinV=-KappaTildeS2BinV;
plot(VBinV,KappaTildeS2BinV,'ro');


