%% Preloading
%Define the physical constant
mli=9.988346*10^-27; %kg
hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant
omega=23.9*2*pi; %in rad/s
pixellength=0.7*10^-6*3; %in m
sigma0=0.215*10^-12/2; %in m^2
%load all the functions
addpath('../../Library');

Fudge=1.458;
kB=1.380e-23;
Nsat=770;
load('/Users/Zhenjie/Data/Processed/2017-05-14/2017-05-14Set3.mat')
warning ('off','all')
%% Get the profile for all of them
nS1list={};
ZS1list={};
Z0S1list=[];
ZS1sortlist={};
EFS1List={};
for i=1:length(imglistS1Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( imglistS1Final{i},'ROI1',[270,265,456,600],...
    'ROI2',[320,388,456,475],'TailRange',[325,580],'ShowOutline',i==1,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',0,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',1,'IfLookUpTable',1,'Zaveraging',0);%'TailRange',[180,380],
    nS1list=[nS1list;nsort/1e18];%
    Z0S1=zcor.z0*pixellength/1e-6;
    Z0S1list=[Z0S1list;Z0S1];
    ZS1=Zsort/1e-6+Z0S1;
    ZS1list=[ZS1list;ZS1];
    ZS1sortlist=[ZS1sortlist;Zsort];
    EFS1List=[EFS1List,EFS1];
end
%% 'ROI2',[300,515,712,930],
nS2list={};
ZS2list={};
Z0S2list=[];
for i=1:length(imglistS2Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( imglistS2Final{i},'ROI1',[270,265,456,600],...
    'ROI2',[270,418,456,450],'TailRange',[325,580],'ShowOutline',i==1,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',0,'OutlineIntrapolate',1,'Zaveraging',0,'IfExtrapolateAngle',1);%
    nS2list=[nS2list;nsort/1e18];
    Z0S2=zcor.z0*pixellength/1e-6;
    Z0S2list=[Z0S2list;Z0S2];
    ZS2=Zsort/1e-6+Z0S2;
    ZS2list=[ZS2list;ZS2];
end
%%
nS1=[];
ZS1=[];
for i=1:length(nS1list)
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
plot(ZS1,nS1,'b.',ZS2,nS2,'r.')
%% Bin the profile with Z
Nbin=250;
ZGrid=linspace(-400,400,Nbin+1);
[ZS1BinZ,nS1Bin]=BinGrid(ZS1,nS1,ZGrid,0);
[ZS2BinZ,nS2Bin]=BinGrid(ZS2,nS2,ZGrid,0);
%nS1Bin=TailTailor(nS1Bin,ZS1BinZ,-170,170);

plot(ZS1BinZ,nS1Bin,'b',ZS2BinZ,nS2Bin,'r')
mask1=isnan(nS1Bin);
nS1Bin(mask1)=[];
ZS1BinZ(mask1)=[];
mask2=isnan(nS2Bin);
ZS2BinZ(mask2)=[];
nS2Bin(mask2)=[];

ZBinZ=ZS1BinZ;
nS1BinZ=nS1Bin;
nS2BinZ=interp1(ZS2BinZ,nS2Bin,ZBinZ,'spline');
VBinZ=0.5*mli*omega^2*(ZBinZ*1e-6).^2/hh;
plot(ZBinZ,nS1BinZ,'b.',ZBinZ,nS2BinZ,'r.');
xlabel('z(um)');ylabel('n (um^{-3})')
%plot(VBinZ,nS1BinZ,'b.',VBinZ,nS2BinZ,'r.');
ntotBinZ=nS1BinZ+nS2BinZ;
EFS1BinZ=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS1BinZ)*1e18).^(2/3).*sign(nS1BinZ)/hh;
EFS2BinZ=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS2BinZ)*1e18).^(2/3).*sign(nS2BinZ)/hh;

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
EFS2=(1/(2*mli))*hbar^2*(6*pi^2*abs(nS2)*1e18).^(2/3).*sign(nS2)/hh;

Nbin=150;
VGrid=linspace(0,1.0e4,Nbin+1);
[VS1BinV,EFS1BinV,VS1ErrBinV,EFS1ErrBinV]=BinGrid(VS1,EFS1,VGrid,0);
mask=isnan(VS1BinV);
VS1BinV(mask)=[];
EFS1BinV(mask)=[];
VS1ErrBinV(mask)=[];
EFS1ErrBinV(mask)=[];

[VS2BinV,EFS2BinV1,VS2ErrBinV1,EFS2ErrBinV1]=BinGrid(VS2,EFS2,VGrid,0);
mask=isnan(VS2BinV);
VS2BinV(mask)=[];
EFS2BinV1(mask)=[];
VS2ErrBinV1(mask)=[];
EFS2ErrBinV1(mask)=[];


Vth=7000;
VtailS1=VS1BinV(VS1BinV>Vth);
EFtailS1=EFS1BinV(VS1BinV>Vth);

Vtail2=VS2BinV(VS2BinV>Vth);
EFtailS2=EFS2BinV1(VS2BinV>Vth);

EFS1BinV=EFS1BinV-mean(EFtailS1);
EFS2BinV1=EFS2BinV1-mean(EFtailS2);


VBinV=VS1BinV;
%nS1BG=mean(nS1BinV(VBinV>8000));
EFS2BinV=interp1(VS2BinV,EFS2BinV1,VBinV);
EFS2ErrBinV=interp1(VS2BinV,EFS2ErrBinV1,VBinV);
nS1BinV=(2*mli*abs(EFS1BinV)*hh/hbar^2).^(3/2)/(6*pi^2)/1e18.*sign(EFS1BinV);
nS2BinV=(2*mli*abs(EFS2BinV)*hh/hbar^2).^(3/2)/(6*pi^2)/1e18.*sign(EFS2BinV);

%nS2BinV=nS2BinV-nS2BG;


ZBinV=sqrt(2*VBinV*hh/(mli*omega^2))/1e-6;



%plot(VBinV,nS1BinV,'r.',VBinV,nS2BinV,'b.');

Vprofile.V=VBinV;
Vprofile.nS1=nS1BinV;
Vprofile.nS2=nS2BinV;
Vprofile.EFS1=EFS1BinV;
Vprofile.EFS2=EFS2BinV;
Vprofile.Z=ZBinV;

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);
plot(VBinV,EFS1BinV,'r.','Displayname','Majority')
hold on
plot(VBinV,EFS2BinV,'b.','Displayname','Minority');
hold off
xlim([0,10000]);
legend('show')
xlabel('U(Hz)');ylabel('E_F(Hz)')
%%
%errorbar(VS1BinV,EFS1BinV,EFS1ErrBinV);
errorbar(VBinV,EFS2BinV,EFS2ErrBinV);
%%
Vth=2100;
mask=VBinV>Vth;

Vfit=VBinV(mask);
nS1fit=nS1BinV(mask);

P=IdealGasFit( Vfit,nS1fit,[5000,500,1] );
P(3)
%%
P=IdealGasFitwoFudge( Vfit,nS1fit,[3000,500] )

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
KappaFitZ(KappaFitZ<0)=0;
% get the kappa from image profile
[KappaTildeS1BinV,KappaTildeS1ErrBinV]=FiniteD( VBinV,VBinV*0,EFS1BinV,EFS1ErrBinV,5);
KappaTildeS1BinV=-KappaTildeS1BinV;

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);
plot(VBinV,KappaTildeS1BinV,'r.');
hold on
plot(VBinZ,KappaFitZ,'k-')
hold off
title('Majority Compressibility with fit')
xlim([0,10e3]);ylim([-0.5,2.2])
ylabel('\kappa/\kappa_0');
xlabel('U (Hz)');
%%

[KappaTildeS2BinV,KappaTildeS2ErrBinV]=FiniteD( VBinV,VBinV*0,EFS2BinV,EFS2ErrBinV,5);
KappaTildeS2BinV=-KappaTildeS2BinV;
plot(VBinV,KappaTildeS2BinV,'ro');
xlim([0,10e3])
ylabel('\kappa/\kappa_0');
xlabel('U (Hz)');
%%
errorbar(VBinV,KappaTildeS2BinV,KappaTildeS2ErrBinV)
xlim([0,10e3])
ylabel('\kappa/\kappa_0');
xlabel('U (Hz)');

%% Polaron EoS
meff=1.25
TTildeS2=T_trap./EFS2BinV;
mask1=TTildeS2>0;
mask2=TTildeS2<2;
mask=mask1 & mask2;
TEoSS2=TTildeS2(mask);
KappaEoSS2=KappaTildeS2BinV(mask)/1.615;
KappaEoSErrS=KappaTildeS2ErrBinV(mask)/1.615;

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,2.7,2.1]);
errorbar(TEoSS2,KappaEoSS2,KappaEoSErrS,'b.','markersize',10,'displayname','Data');
hold on 
plot(TTildeT,KappaTildeT,'k-','displayname','m*/m=1');
plot(TTildeT/meff,KappaTildeT*meff,'k--','displayname',['m*/m=',num2str(meff)]);
hold off
xlim([0,2]);%ylim([0.2,1.5])
ylabel('\kappa/\kappa_0');
xlabel('T/T_F');
legend('show')

%%
n_ideal_BinV=IdealGasnV( [mu_trap_S1,T_trap],VBinV );

plot(VBinV,n_ideal_BinV,'k-','linewidth',1,'DisplayName','ideal');
hold on
%plot(VBinZ,nS1BinZ,'bo','DisplayName','Majority');
plot(VBinV,nS1BinV,'r.','DisplayName','Majority');
plot(VBinV,nS2BinV,'b.','DisplayName','Minority');
plot(VBinV,n_ideal_BinV+0.6150*nS2BinV,'g-','linewidth',1,'DisplayName','ENS EoS, A=-0.615');
hold off
xlim([0,8000])
legend show

ylabel('n (um^{-3})');
xlabel('U (Hz)');


%%
%Imbalance=max(nS2BinV)/max(nS1BinV)
Imbalance=(max(nS1BinV)-max(nS2BinV))/(max(nS1BinV)+max(nS2BinV))
%% Create the unfolded kappa vs Z

ZList=[];
KappaList=[];
Iteration=10;
for i=1:Iteration*length(ZS1sortlist)
    % first bin the data
    %Zgrid=linspace(-200,200,Nbin+1);
    %Zgrid=linspace(-sqrt(200),sqrt(200),Nbin+1);Zgrid=sign(Zgrid).*Zgrid.^2;
    j=mod(i,length(ZS1sortlist));
    if j==0
        j=length(ZS1sortlist);
    end
    k=round((i-j)/length(ZS1sortlist))+1;
    Nbin=600;
    Z_vec=ZS1sortlist{j};
    Tgrid0=linspace(sign(min(Z_vec))*min(Z_vec).^2,sign(max(Z_vec))*max(Z_vec).^2,Nbin+1);
    DeltaT=abs(Tgrid0(2)-Tgrid0(1));
    Toffset=linspace(-0.5*DeltaT,0.5*DeltaT,Iteration);
    Tgrid=Tgrid0+Toffset(k);
    Zgrid=sign(Tgrid).*sqrt(abs(Tgrid));
    
    [ZBintemp,EFBintemp,~,~]=BinGrid(ZS1sortlist{j},EFS1List{j}/hh,Zgrid,0);

    % then divide the data into two group Z>0 and Z<0
    Ztempplus=ZBintemp(ZBintemp>=0);EFtempplus=EFBintemp(ZBintemp>=0);
    Ztempminus=ZBintemp(ZBintemp<0);EFtempminus=EFBintemp(ZBintemp<0);
    Vtempplus=0.5*mli*omega^2*(Ztempplus).^2/hh;
    Vtempminus=0.5*mli*omega^2*(Ztempminus).^2/hh;
    [kappatempplus,~] = FiniteD( Vtempplus,Vtempplus*0,EFtempplus,EFtempplus*0,2);
    kappatempplus=-kappatempplus;
    [kappatempminus,~] = FiniteD( Vtempminus,Vtempminus*0,EFtempminus,EFtempminus*0,2);
    kappatempminus=-kappatempminus;
 
    pointsskip=0;
    
    Ztemp=[Ztempminus(1:end-pointsskip),Ztempplus(pointsskip+1:end)];
    kappatemp=[kappatempminus(1:end-pointsskip),kappatempplus(pointsskip+1:end)];
    ZList=[ZList,Ztemp];
    KappaList=[KappaList,kappatemp];
    disp(i)
end
ZList=ZList/1e-6;
%%
plot(ZList,KappaList,'r.','markersize',20)
ylim([-1,2]);

%%
Nbin=100;
Zgrid1=linspace(-(250^2),250^2,101+1);
Zgrid1=sqrt(abs(Zgrid1)).*sign(Zgrid1);
Zgrid2=linspace(-350,350,Nbin+1);

Zgrid=sort([Zgrid1(abs(Zgrid1)<100),Zgrid2(abs(Zgrid2)>100)]);
%Zgrid=Zgrid2;

[ZBinK,KappaBinK,~,~,~,KappaBinErrK]=BinGrid(ZList,KappaList,Zgrid,2);
KappaBinErrK=KappaBinErrK/sqrt(length(ZS1sortlist));
errorbar(ZBinK,KappaBinK,KappaBinErrK,'r.','markersize',20);
ylim([-0.2,3.5]);
%% Kappa vs Z plot

Vth=1700;
Zth=sqrt(2*Vth*hh/(mli*omega^2))/1e-6;

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
plot( ZS1BinZ,KappaFitZ,'Color','k')
hold on
%plot(ZBin,KappaBin,'.','markersize',20,'color',[49,115,255]/255,'linewidth',1)
errorbar1derr_Z( ZBinK,KappaBinK,KappaBinErrK,'MarkerEdgeColor',[49,115,255]/255,'ErrLineWidth',0.5,'ErrBarColor',[36,85,189]/255,'LineStyle','none','Markersize',5)

line([Zth,Zth],[-20,20],'linewidth',0.5,'color','k')
line([-Zth,-Zth],[-20,20],'linewidth',0.5,'color','k')

hold off
ylim([-0.2,3.8]);xlim([-250,250]);
xlabel('z(um)');ylabel('KappaTilde')
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'box','on','Xtick',[-200,-100,0,100,200]);

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.5,0.6]);
mask=abs(ZBinZ)<170;
%plot(ZS1grid/1e-6,TTilde,'-')
plot(ZBinZ(mask),T_trap./EFS1BinZ(mask),'-','color',[36,85,189]/255,'linewidth',1)
hold on
line([Zth,Zth],[-100,100],'color','k','linewidth',0.5)
line([-Zth,-Zth],[-100,100],'color','k','linewidth',0.5)

hold off
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0 0.1 0.2 0.4 0.6],'Xtick',[-200,-100,0,100,200])
ylabel('T/T_F');xlabel('Z');
xlim([-250,250]);ylim([0,0.3]);

%%

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);
plot(ZBinZ,nS2BinZ,'color',[201,67,52]/255,'linewidth',1);
hold on
plot(ZBinZ,nS1BinZ,'color',[36,85,189]/255,'linewidth',1)
line([Zth,Zth],[-1000,1000],'linewidth',0.5,'color','k');
line([-Zth,-Zth],[-1000,1000],'linewidth',0.5,'color','k');
%plot(ZS1BinZ,nS2BinZ./nS1BinZ,'g-','linewidth',1)
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0,0.1,0.2,0.3],'Xtick',[-200,-100,0,100,200])
hold off

xlim([-250,250]);xlabel('z(um)')
ylim([-0.05,0.4]);ylabel('n(um^{-3})')

%%
ROI=[280,155,407,462];
AverageCountsA=mean(mean(imglistS1Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),2)-imglistS1Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),3)))
AverageCountsB=mean(mean(imglistS2Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),2)-imglistS2Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),3)))
