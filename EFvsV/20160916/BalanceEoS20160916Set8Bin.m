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
Fudge=1.9442;
kB=1.380e-23;
Nsat=770;
load('/Users/Zhenjie/Data/Processed/2016-09-16/2016-09-16Set8Bin.mat')
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
    'ROI2',[270,400,456,450],'TailRange',[325,580],'ShowOutline',i==1,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',0,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',1,'IfLookUpTable',1);
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
ZS2sortlist={};
EFS2List={};
for i=1:length(imglistS1Final)
    disp(i);
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS2,P,zcor,Vsel]=EOS_Online( imglistS2Final{i},'ROI1',[270,265,456,600],...
    'ROI2',[270,400,456,450],'TailRange',[325,580],'ShowOutline',i==1,'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',Fudge,'BGSubtraction',0,'IfFitExpTail',0,'Nsat',Nsat,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',pixellength,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',1,'IfLookUpTable',1);
    nS2list=[nS2list;nsort/1e18];
    Z0S2=zcor.z0*pixellength/1e-6;
    Z0S2list=[Z0S2list;Z0S2];
    ZS2=Zsort/1e-6+Z0S2;
    ZS2list=[ZS2list;ZS2];
    ZS2sortlist=[ZS2sortlist;Zsort];
    EFS2List=[EFS2List,EFS2];
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
for i=1:length(nS2list)
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
ZGrid=linspace(-300,300,Nbin+1);
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

Nbin=200;
VGrid=linspace(0,0.8e4,Nbin+1);
[VS1BinV,EFS1BinV]=BinGrid(VS1,EFS1,VGrid,0);
mask=isnan(VS1BinV);
VS1BinV(mask)=[];
EFS1BinV(mask)=[];

[VS2BinV,EFS2BinV1]=BinGrid(VS2,EFS2,VGrid,0);
mask=isnan(VS2BinV);
VS2BinV(mask)=[];
EFS2BinV1(mask)=[];

Vth=3000;
VtailS1=VS1BinV(VS1BinV>Vth);
EFtailS1=EFS1BinV(VS1BinV>Vth);

Vtail2=VS2BinV(VS2BinV>Vth);
EFtailS2=EFS2BinV1(VS2BinV>Vth);

EFS1BinV=EFS1BinV-mean(EFtailS1);
EFS2BinV1=EFS2BinV1-mean(EFtailS2);


VBinV=VS1BinV;
%nS1BG=mean(nS1BinV(VBinV>8000));
EFS2BinV=interp1(VS2BinV,EFS2BinV1,VBinV);
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
xlim([0,8000]);
legend('show')
xlabel('U(Hz)');ylabel('E_F(Hz)')
%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);
% get the kappa from image profile
[KappaTildeS1BinV,~]=FiniteD( VBinV,VBinV*0,EFS1BinV,EFS1BinV*0,5);
KappaTildeS1BinV=-KappaTildeS1BinV;
plot(VBinV,KappaTildeS1BinV,'b.');
title('Compressibility S1')
xlim([0,7e3]);ylim([-0.5,3.5])
ylabel('\kappa/\kappa_0');
xlabel('U (Hz)');
%%

[KappaTildeS2BinV,~]=FiniteD( VBinV,VBinV*0,EFS2BinV,EFS2BinV*0,4);
KappaTildeS2BinV=-KappaTildeS2BinV;
plot(VBinV,KappaTildeS2BinV,'ro');
xlim([0,10e3])
ylabel('\kappa/\kappa_0');
xlabel('U (Hz)');
%%
hold on
plot(VBinV,nS1BinV,'b.','DisplayName','Majority');
plot(VBinV,nS2BinV,'r.','DisplayName','Minority');
hold off
xlim([0,6000])
legend show

ylabel('n (um^{-3})');
xlabel('U (Hz)');
%%
Imbalance=max(nS2BinV)/max(nS1BinV)
%% Get the reduced pressure
np=real(EFS1BinV.^(3/2));
P=np;
VthP=3500;
np(VBinV>VthP)=0;

for i=1:(length(P)-1)
    P(i)=real(trapz(VBinV(i:end),np(i:end)));
end
P0=0.4*np.*EFS1BinV;
P1T=P./P0;

plot(VBinV,P1T);ylim([0,4])
%% fit the pressure with Mark's EoS to get temperature
load('/Users/Zhenjie/Data/Processed/Mark/MarkEoS.mat')
Vth1=500;
Vth2=1400;
mask1=VBinV>Vth1;mask2=VBinV<Vth2;
mask=mask1 & mask2;
VSample=VBinV(mask);
EFSample=EFS1BinV(mask);
PSample=P1T(mask);
TtildeSample=interp1(PTildeMark,TTildeMark,PSample);
TSample=TtildeSample.*EFSample;
scatter(VSample,TSample);
ylim([0,max(TSample)]);
T=mean(TSample)
plot(VSample,TSample,'r.')
ylabel('k_B T (Hz)');xlabel('V(Hz)');
title('Temperature get from P/P_0');
ylim([0,1.2*max(TSample)])
%%
%% Get the mu for the trap from mark's EoS

T_trap=T;
Vmax=1500;
mask=VS1BinV<Vmax;
Vfitmu=VS1BinV(mask);
EFfitmu=EFS1BinV(mask);

TTildefitmu=T_trap./EFfitmu;
mu_EF_U_fitmu=interp1(TTildeMark,mu_EF_Mark,TTildefitmu,'spline');

mu_local_fitmu=mu_EF_U_fitmu.*EFfitmu;
mu_fitmu=mu_local_fitmu+Vfitmu;

mask1=TTildefitmu<min(TTildeMark);
mask2=TTildefitmu>max(TTildeMark);
mask=mask1 | mask2;

Vfitmu(mask)=[];
mu_fitmu(mask)=[];

scatter(Vfitmu,mu_fitmu);

mu_trap=mean(mu_fitmu)


%% Create the unfolded kappa vs Z

ZList=[];
KappaList=[];
Iteration=10;
Zoffset=0e-6;

for i=1:Iteration*length(ZS1sortlist)
    % first bin the data
    %Zgrid=linspace(-200,200,Nbin+1);
    %Zgrid=linspace(-sqrt(200),sqrt(200),Nbin+1);Zgrid=sign(Zgrid).*Zgrid.^2;
    
    j=mod(i,length(ZS1sortlist));
    if j==0
        j=length(ZS1sortlist);
    end
    
    k=round((i-j)/length(ZS1sortlist))+1;
    Nbin=2500;
    Z_vec=ZS1sortlist{j}-Zoffset;
    Tgrid0=linspace(sign(min(Z_vec))*min(Z_vec).^2,sign(max(Z_vec))*max(Z_vec).^2,Nbin+1);
    DeltaT=abs(Tgrid0(2)-Tgrid0(1));
    Toffset=linspace(-0.5*DeltaT,0.5*DeltaT,Iteration);
    Tgrid=Tgrid0+Toffset(k);
    Zgrid=sign(Tgrid).*sqrt(abs(Tgrid));
    
    [ZBintemp,EFBintemp,~,~]=BinGrid(Z_vec,EFS1List{j}/hh,Zgrid,0);

    % then divide the data into two group Z>0 and Z<0
    Ztempplus=ZBintemp(ZBintemp>=0);EFtempplus=EFBintemp(ZBintemp>=0);
    Ztempminus=ZBintemp(ZBintemp<0);EFtempminus=EFBintemp(ZBintemp<0);
    Vtempplus=0.5*mli*omega^2*(Ztempplus).^2/hh;
    Vtempminus=0.5*mli*omega^2*(Ztempminus).^2/hh;
    [kappatempplus,~] = FiniteD( Vtempplus,Vtempplus*0,EFtempplus,EFtempplus*0,1);
    kappatempplus=-kappatempplus;
    [kappatempminus,~] = FiniteD( Vtempminus,Vtempminus*0,EFtempminus,EFtempminus*0,1);
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
ylim([-1,3.5]);

%%
Nbin=80;
Zgrid1=linspace(-(250^2),250^2,140);
Zgrid1=sqrt(abs(Zgrid1)).*sign(Zgrid1);
Zgrid2=linspace(-250,250,Nbin+1);

Zgrid=sort([Zgrid1(abs(Zgrid1)<120),Zgrid2(abs(Zgrid2)>120)]);
%Zgrid=Zgrid2;

[ZBinK,KappaBinK,~,~,~,KappaBinErrK]=BinGrid(ZList,KappaList,Zgrid,2);
KappaBinErrK=KappaBinErrK/sqrt(length(ZS1sortlist));
errorbar(ZBinK,KappaBinK,KappaBinErrK,'r.','markersize',20);
ylim([-0.2,3.5]);


%% get a fit for kappa/kappa_0 V2
%Zfit=linspace(-200,200,100);
Vgrid=0.5*mli*omega^2*(linspace(0,200,100)*1e-6).^2/hh;

mu_local=mu_trap-Vgrid;
Beta_mu_local=mu_local/T_trap;

% KappaFit=interp1(mu_T_Mark,KappaTildeMark,Beta_mu_local,'linear');
Vfit=interp1(Beta_mu_local,Vgrid,mu_T_Mark,'linear');
Zfit=sqrt(Vfit*hh*2/(mli*omega^2))/1e-6;
Zfit=[-Zfit;Zfit];
KappaFit=[KappaTildeMark;KappaTildeMark];
KappaFitErr=[dKappaTildeMark;dKappaTildeMark];
% mask1=Beta_mu_local<min(mu_T_Mark);
% mask2=Beta_mu_local>max(mu_T_Mark);
% KappaFit(mask1)=0;
% KappaFit(mask2)=KappaTildeMark(1);

% KappaFit(mask1)=nan;
% KappaFit(mask2)=nan;
% KappaFitErr(mask1)=nan;
% KappaFitErr(mask2)=nan;

plot(Zfit,KappaFit,'b.');

ylim([0,5]);

%% Get the critical point
mask=abs(ZBinZ)<95;
ZTshow=ZBinZ(mask);
Tshow=T_trap./EFS1BinZ(mask);

Zcr1=interp1(Tshow(ZTshow>0),ZTshow(ZTshow>0),0.167);
Zcr2=interp1(Tshow(ZTshow<0),ZTshow(ZTshow<0),0.167);

%% plot n vs z

figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);

hold on
plot(ZBinZ,nS2BinZ,'color',[201,67,52]/255,'linewidth',1);
plot(ZBinZ,nS1BinZ,'color',[36,85,189]/255,'linewidth',1);
% line([Zcr1,Zcr1],[-10,10],'linewidth',1,'color',[201,67,52]/255);
% line([Zcr2,Zcr2],[-10,10],'linewidth',1,'color',[201,67,52]/255);
line([Zcr1,Zcr1],[-20,20],'linewidth',0.5,'color','k')
line([Zcr2,Zcr2],[-20,20],'linewidth',0.5,'color','k')
box on
hold off
ylim([-0.05,0.2]);xlim([-250,250])
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0 0.1 0.2 0.3],'Xtick',[-200,-100,0,100,200])
xlabel('Z (um)');ylabel('n(um^{-3})')

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.5,0.6]);
plot(ZBinZ,nS2BinZ,'color',[201,67,52]/255,'linewidth',1);
hold on
plot(ZBinZ,nS1BinZ,'color',[36,85,189]/255,'linewidth',1)
line([Zcr1,Zcr1],[-1000,1000],'linewidth',0.5,'color','k');
line([Zcr2,Zcr2],[-1000,1000],'linewidth',0.5,'color','k');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0,0.1,0.2],'Xtick',[-200,-100,0,100,200])
hold off

xlim([-250,250]);
ylim([-0.05,0.4]);
%% T/TF plot
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.5,0.6]);
%plot(ZS1grid/1e-6,TTilde,'-')
plot(ZTshow,Tshow,'-','color',[36,85,189]/255,'linewidth',1)
line([Zcr1,Zcr1],[-20,20],'linewidth',0.5,'color','k')
line([Zcr2,Zcr2],[-20,20],'linewidth',0.5,'color','k')
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0 0.1 0.2],'Xtick',[-200,-100,0,100,200])
ylabel('T/T_F');xlabel('Z');
xlim([-250,250]);ylim([0,0.3]);

%% Kappa vs Z plot
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);%plot(Zfit,KappaFit,'linewidth',1,'color','k');

%errorbar1derr_Z( Zfit,KappaFit,KappaFitErr,'MarkerEdgeColor',[65,64,66]/255,'ErrLineWidth',0.5,'ErrBarColor','k','LineStyle','none','Markersize',3 )
%shadedErrorBar(Zfit(Zfit<0),KappaFit(Zfit<0),KappaFitErr(Zfit<0))
hold on
%shadedErrorBar(Zfit(Zfit>0),KappaFit(Zfit>0),KappaFitErr(Zfit>0))
%plot(Z_kShow,KappaUShow,'.','markersize',5,'color',[49,115,255]/255,'linewidth',1)
errorbar1derr_Z( ZBinK,KappaBinK,KappaBinErrK,'MarkerEdgeColor',[49,115,255]/255,'ErrLineWidth',0.5,'ErrBarColor',[36,85,189]/255,'LineStyle','none','Markersize',5 )
line([Zcr1,Zcr1],[-20,20],'linewidth',0.5,'color','k')
line([Zcr2,Zcr2],[-20,20],'linewidth',0.5,'color','k')
line([-250,250],[1/0.376,1/0.376],'linewidth',0.5,'color','k');
line([-250,250],[1/0.42,1/0.42],'linewidth',0.5,'color','k');
hold off
ylim([-0.2,3.8]);xlim([-250,250])
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0],'Ytick',[0 1 2 3],'Xtick',[-200,-100,0,100,200],'box','on')
ylabel('\kappa/\kappa_0');xlabel('z (um)')

%%
ROI=[280,155,407,462];
AverageCountsA=mean(mean(imglistS1Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),2)-imglistS1Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),3)))
AverageCountsB=mean(mean(imglistS2Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),2)-imglistS2Final{1}(ROI(2):ROI(4),ROI(1):ROI(3),3)))
