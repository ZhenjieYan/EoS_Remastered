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
Fudge=1.1935;
%Fudge=1.7087;
kB=1.380e-23;
Nsat=770;
[ KappaTildeT, PTildeT, TTildeT, CV_NkBT , beta_mu_vecT ,Z_vecT ] = IdealFermiEOS();
load('/Users/Zhenjie/Data/Processed/2016-10-20/2016-10-20Set10Bin.mat')
warning ('off','all')
% imglistS1Final=imglistS1Final(1:4);
% imglistS2Final=imglistS2Final(1:4);
%% Get the profile for all of them
nS1list={};
ZS1list={};
Z0S1list=[];
ZS1sortlist={};
EFS1List={};
for i=1:length(imglistS1Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel]=EOS_Online( imglistS1Final{i},'ROI1',[220,265,426,630],...
    'ROI2',[220,400,426,490],'TailRange',[325,580],'ShowOutline',i==1,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',0,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',1,'IfLookUpTable',1,'Zaveraging',0,'BoxShape','Square');%'TailRange',[180,380],
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
EFS2List={};
for i=1:length(imglistS2Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS2,P,zcor,Vsel]=EOS_Online( imglistS2Final{i},'ROI1',[220,265,426,630],...
    'ROI2',[270,430,456,450],'TailRange',[325,580],'ShowOutline',i==1,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',1,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',1,'OutlineIntrapolate',1,'Zaveraging',0,'IfExtrapolateAngle',1,'BoxShape','Square');%
    nS2list=[nS2list;nsort/1e18];
    Z0S2=zcor.z0*pixellength/1e-6;
    Z0S2list=[Z0S2list;Z0S2];
    ZS2=Zsort/1e-6+Z0S2;
    ZS2list=[ZS2list;ZS2];
    EFS2List=[EFS2List,EFS2];
end
%%
nS1=[];
ZS1=[];
EFS1=[];
for i=1:length(nS1list)
    nS1=[nS1;nS1list{i}];
    ZS1=[ZS1;ZS1list{i}-mean(Z0S1list)];
    EFS1=[EFS1;EFS1List{i}/hh];
end
nS2=[];
ZS2=[];
EFS2=[];
for i=1:length(imglistS2Final)
    nS2=[nS2;nS2list{i}];
    ZS2=[ZS2;ZS2list{i}-mean(Z0S2list)];
    EFS2=[EFS2;EFS2List{i}/hh];
end

UnbinnedProfile.nS1=nS1;
UnbinnedProfile.nS2=nS2;
UnbinnedProfile.ZS1=ZS1;
UnbinnedProfile.ZS2=ZS2;
plot(ZS1,EFS1,'b.',ZS2,EFS2,'r.')
%% Bin the profile with Z for EF
Nbin=250;
ZGrid=linspace(-400,400,Nbin+1);
[ZS1BinEF,EFS1BinEF]=BinGrid(ZS1,EFS1,ZGrid,0);
[ZS2BinEF,EFS2BinEF]=BinGrid(ZS2,EFS2,ZGrid,0);
plot(ZS1BinEF,EFS1BinEF,'b.',ZS2BinEF,EFS2BinEF,'r.')
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
[VS1BinV,EFS1BinV]=BinGrid(VS1,EFS1,VGrid,0);
mask=isnan(VS1BinV);
VS1BinV(mask)=[];
EFS1BinV(mask)=[];
[VS2BinV,EFS2BinV1,VS2ErrBinV1,EFS2ErrBinV1]=BinGrid(VS2,EFS2,VGrid,0);
mask=isnan(VS2BinV);
VS2BinV(mask)=[];
EFS2BinV1(mask)=[];
VS2ErrBinV1(mask)=[];
EFS2ErrBinV1(mask)=[];

Vth=4000;
VtailS1=VS1BinV(VS1BinV>Vth);
EFtailS1=EFS1BinV(VS1BinV>Vth);

Vtail2=VS2BinV(VS2BinV>Vth);
EFtailS2=EFS2BinV1(VS2BinV>Vth);

EFS1BinV=EFS1BinV-mean(EFtailS1);
EFS2BinV1=EFS2BinV1-mean(EFtailS2);



Vth=7000;
VtailS1=VS1BinV(VS1BinV>Vth);
EFtailS1=EFS1BinV(VS1BinV>Vth);

Vtail2=VS2BinV(VS2BinV>Vth);
EFtailS2=EFS2BinV1(VS2BinV>Vth);

EFS1BinV=EFS1BinV-mean(EFtailS1);
EFS2BinV1=EFS2BinV1-mean(EFtailS2);


VBinV=VS1BinV;
%nS1BG=mean(nS1BinV(VBinV>8000));
EFS2BinV=interp1(VS2BinV,EFS2BinV1,VBinV,'cube');
EFS2ErrBinV=interp1(VS2BinV,EFS2ErrBinV1,VBinV,'cube');
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
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,2.4]);
plot(VBinV,EFS1BinV,'b.','Displayname','Majority')
hold on
plot(VBinV,EFS2BinV,'r.','Displayname','Minority');
hold off
xlim([0,10000]);
legend('show')
xlabel('U(Hz)');ylabel('E_F(Hz)')
%%
Vth=1200;
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
[KappaTildeS1BinV,~]=FiniteD( VBinV,VBinV*0,EFS1BinV,EFS1BinV*0,9);
KappaTildeS1BinV=-KappaTildeS1BinV;

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,2.4]);
plot(VBinV,KappaTildeS1BinV,'b.');
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
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);

n_ideal_BinV=IdealGasnV( [mu_trap_S1,T_trap],VBinV );
plot(VBinV,n_ideal_BinV,'k-','linewidth',1,'DisplayName','ideal');
hold on
plot(VBinV,nS1BinV,'b.','DisplayName','Majority');
plot(VBinV,nS2BinV,'r.','DisplayName','Minority');
plot(VBinV,n_ideal_BinV+0.6150*nS2BinV,'g-','linewidth',1,'DisplayName','ENS EoS, A=-0.615');
hold off
xlim([0,10000])
legend show

ylabel('n (um^{-3})');
xlabel('U (Hz)');
%% Polaron EoS
meff=1.2
TTildeS2=T_trap./EFS2BinV;
mask1=TTildeS2>0;
mask2=TTildeS2<2;
mask=mask1 & mask2;
TEoSS2=TTildeS2(mask);
KappaEoSS2=KappaTildeS2BinV(mask)/1.615;
KappaEoSErrS=KappaTildeS2ErrBinV(mask)/1.615;



figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
errorbar(TEoSS2,KappaEoSS2,KappaEoSErrS,'b.','markersize',10,'displayname','Data');
hold on 
plot(TTildeT,KappaTildeT,'k-','displayname','m*/m=1');
plot(TTildeT/meff,KappaTildeT*meff,'k--','displayname',['m*/m=',num2str(meff)]);
hold off
xlim([0,2]);ylim([0.2,1.5])
ylabel('\kappa/\kappa_0');
xlabel('T/T_F');
legend('show')


%%
Imbalance=max(nS2BinV)/max(nS1BinV)
%%
P=TFfit( nS2BinZ,ZBinZ,100);
ZTF=P(3);

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
    Nbin=1200;
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
    [kappatempplus,~] = FiniteD( Vtempplus,Vtempplus*0,EFtempplus,EFtempplus*0,3);
    kappatempplus=-kappatempplus;
    [kappatempminus,~] = FiniteD( Vtempminus,Vtempminus*0,EFtempminus,EFtempminus*0,3);
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

line([ZTF,ZTF],[-20,20],'linewidth',0.5,'color','k')
line([-ZTF,-ZTF],[-20,20],'linewidth',0.5,'color','k')

hold off
ylim([-0.2,3.8]);xlim([-250,250]);
xlabel('z(um)');ylabel('KappaTilde')
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'box','on','Xtick',[-200,-100,0,100,200]);

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,0.63]);
mask=abs(ZBinZ)<170;
%plot(ZS1grid/1e-6,TTilde,'-')
plot(ZBinZ(mask),T_trap./EFS1BinZ(mask),'-','color',[36,85,189]/255,'linewidth',1)
hold on
% line([Zth,Zth],[-100,100],'color','k','linewidth',0.5)
% line([-Zth,-Zth],[-100,100],'color','k','linewidth',0.5)

line([ZTF,ZTF],[-100,100],'color','k','linewidth',0.5)
line([-ZTF,-ZTF],[-100,100],'color','k','linewidth',0.5)

hold off
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0 0.1 0.2 0.3 0.4 ],'Xtick',[-200,-100,0,100,200])
ylabel('T/T_F');xlabel('Z');
xlim([-250,250]);ylim([0,0.4]);

%%


figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
plot(ZBinZ,nS2BinZ,'color',[201,67,52]/255,'linewidth',1);
hold on
plot(ZBinZ,nS1BinZ,'color',[36,85,189]/255,'linewidth',1)
% line([Zth,Zth],[-2000,5000],'linewidth',0.5,'color','k');
% line([-Zth,-Zth],[-2000,5000],'linewidth',0.5,'color','k');


line([ZTF,ZTF],[-100,100],'color','k','linewidth',0.5)
line([-ZTF,-ZTF],[-100,100],'color','k','linewidth',0.5)


set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0,0.1,0.2,0.3],'Xtick',[-200,-100,0,100,200])
hold off

xlim([-250,250]);xlabel('z(um)')
ylim([-0.03,0.27]);ylabel('n(um^{-3})')

%%

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,0.63]);
plot(ZS1BinEF,EFS1BinEF/1e3,'color',[201,67,52]/255,'linewidth',1);
hold on
plot(ZS2BinEF,EFS2BinEF/1e3,'color',[36,85,189]/255,'linewidth',1)
% line([Zth,Zth],[-2000,5000],'linewidth',0.5,'color','k');
% line([-Zth,-Zth],[-2000,5000],'linewidth',0.5,'color','k');


line([ZTF,ZTF],[-100,100],'color','k','linewidth',0.5)
line([-ZTF,-ZTF],[-100,100],'color','k','linewidth',0.5)

set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0,2.5 5],'Xtick',[-200,-100,0,100,200])
hold off

xlim([-250,250]);xlabel('z(um)')
ylim([-0.5,6.5]);ylabel('n(um^{-3})')


%%
ROI=[280,155,407,462];
AverageCountsA=mean(mean(imglistS1Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),2)-imglistS1Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),3)))
AverageCountsB=mean(mean(imglistS2Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),2)-imglistS2Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),3)))
%% Also get the T from P/P0
np=real(EFS1BinV.^(3/2));

VthP=6500;
np(VBinV>VthP)=0;

P=np*0;
for i=1:(length(P)-1)
    P(i)=real(trapz(VBinV(i:end),np(i:end)));
end
P0=0.4*np.*EFS1BinV;
P1T=P./P0;
% scatter(VS1Bin,P1T);
% ylim([0,5]);

Vth1=1300;
Vth2=3000;
mask1=VBinV>Vth1;mask2=VBinV<Vth2;
mask=mask1 & mask2;
VSample=VS1BinV(mask);
EFSample=EFS1BinV(mask);
PSample=P1T(mask);
TtildeSample=interp1(PTildeT,TTildeT,PSample,'spline');
TSample=TtildeSample.*EFSample;
scatter(VSample,TSample);
ylim([0,1000]);
T=mean(TSample)
Tstd=std(TSample)
ylabel('k_B T (Hz)');xlabel('V(Hz)');
title('Temperature get from P/P_0');
T_trap2=T;
%%
TK=T*hh/kB*1e9
TstdK=Tstd*hh/kB*1e9
%% find mu
Vth1=1800;
Vth2=2800;
mask1=VBinV>Vth1;mask2=VBinV<Vth2;
mask=mask1 & mask2;
VSample=VS1BinV(mask);
EFSample=EFS1BinV(mask);
TTildeSample=T_trap2./EFSample;
beta_mu_T_Sample=interp1(TTildeT,beta_mu_T,TTildeSample);
mu_Sample=beta_mu_T_Sample*T_trap2+VSample;
scatter(VSample,mu_Sample);
mu_trap_S1_2=mean(mu_Sample)
%% Get the EoS of ideal fermi gas with the temperature and mu got from P/P_0
mu_EF_T=beta_mu_T.*TTildeT;
mu_localZ=mu_trap_S1_2-VBinZ;
% get the kappa from the known T and mu
beta_mu_localZ2=mu_localZ/T_trap2;
KappaFitZ2=interp1(beta_mu_T,KappaTildeT,beta_mu_localZ2,'spline');
KappaFitZ2(KappaFitZ2<0)=0;
%%
plot(ZS1BinZ,KappaFitZ2)
hold on
plot(ZS1BinZ,KappaFitZ,'r')
hold off

%% Kappa vs Z plot

Vth=1700;
Zth=sqrt(2*Vth*hh/(mli*omega^2))/1e-6;

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,0.84]);
plot( ZS1BinZ,KappaFitZ2,'Color','k')
hold on
% plot(ZS1BinZ,KappaFitZ,'k')
%plot(ZBin,KappaBin,'.','markersize',20,'color',[49,115,255]/255,'linewidth',1)
errorbar1derr_Z( ZBinK,KappaBinK,KappaBinErrK,'MarkerEdgeColor',[49,115,255]/255,'ErrLineWidth',0.5,'ErrBarColor',[36,85,189]/255,'LineStyle','none','Markersize',5)

line([ZTF,ZTF],[-20,20],'linewidth',0.5,'color','k')
line([-ZTF,-ZTF],[-20,20],'linewidth',0.5,'color','k')

hold off
ylim([-0.2,3.8]);xlim([-250,250]);
xlabel('z(um)');ylabel('KappaTilde')
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'box','on','Xtick',[-200,-100,0,100,200]);


%%
% T=15.9 nK
% T_Err=2.0 nK

%%

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,0.63]);
mask=abs(ZBinZ)<170;
%plot(ZS1grid/1e-6,TTilde,'-')
%plot(ZBinZ(mask),T./EFS1BinZ(mask),'-','color',[36,85,189]/255,'linewidth',1)
plot(ZBinZ(mask),(T+Tstd)./EFS1BinZ(mask),'-','color',[36,85,189]/255,'linewidth',1)
hold on
plot(ZBinZ(mask),(T-Tstd)./EFS1BinZ(mask),'-','color',[36,85,189]/255,'linewidth',1)
fill([ZBinZ(mask),fliplr(ZBinZ(mask))],[(T+Tstd)./EFS1BinZ(mask),fliplr((T-Tstd)./EFS1BinZ(mask))],[36,85,189]/255)
% line([Zth,Zth],[-100,100],'color','k','linewidth',0.5)
% line([-Zth,-Zth],[-100,100],'color','k','linewidth',0.5)

line([ZTF,ZTF],[-100,100],'color','k','linewidth',0.5)
line([-ZTF,-ZTF],[-100,100],'color','k','linewidth',0.5)

hold off
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0 0.1 0.2 0.3 0.4 ],'Xtick',[-200,-100,0,100,200])
ylabel('T/T_F');xlabel('Z');
xlim([-250,250]);ylim([0,0.4]);