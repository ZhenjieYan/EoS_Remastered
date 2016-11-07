function [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFsort,P,zcor,Vsel,output] = EOS_Online( Input,varargin )
%The hope is to use this function as an easy way to do online fitting for
%EoS data
%FileName: the name of the file
%Pt:P/P0
%Kt:kappa/kappa_0
%ns:numberdensity
%Vs:potential
%Zs:The z position
IfFourierFilter=0;
CutOffFactor=0.5;
CropROI1=1;
CropROI2=1;
CropTail=1;
ShowOutline=1;
PolyOrder=2;
Points=20;
VrangeFactor=10;
KappaMode=0;
BGSubtraction=0;
ifsmooth=0;
kmin=0.15;
kmax=1.1;
CutOff=2;
IfHalf=0;
FinalPlot=1;
fudge=2.62;
ShowTailor=0;
SelectByPortion=0;
Portion=0.05;
IfTailTailor=1;
Nsat=330;
IfFitExpTail=0;
ExpTailPortion=0.1;
SM=1;
IfBin=0;
BinGridSize=100;
IfSuperSampling=0;
Nmethod=0;
pixellength=1.44*10^-6; %in m
IfCleanImage=0;
OutlineExtrapolate=1;
IfLookUpTable=0;
ImagingTime=10;
Zaveraging=1;
omega=23.9*2*pi; %in rad/s
OutlineIntrapolate=1;
IfExtrapolateAngle=0;
IfReadImg=0;
BoxShape='Circular';
for i =1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            case 'ROI1'
                ROI1=varargin{i+1};
                CropROI1=0;
            case 'ROI2'
                ROI2=varargin{i+1};
                CropROI2=0;
            case 'TailRange'
                Zrange=varargin{i+1};
                zmin=Zrange(1);
                zmax=Zrange(2);
                CropTail=0;
            case 'ShowPlot'
                FinalPlot=varargin{i+1};
            case 'PolynomialPoints'
                Points=varargin{i+1};
            case 'VrangeFactor'
                VrangeFactor=varargin{i+1};
            case 'KappaMode'
                KappaMode=varargin{i+1};
            case 'BGSubtraction'
                BGSubtraction=varargin{i+1};
            case 'ShowOutline'
                ShowOutline=varargin{i+1};
            case 'kmin'
                kmin=varargin{i+1};
            case 'kmax'
                kmax=varargin{i+1};
            case 'smooth'
                ifsmooth=varargin{i+1};
            case 'PolyOrder'
                PolyOrder=varargin{i+1};
            case 'IfHalf'
                IfHalf=varargin{i+1};
            case 'Fudge'
                fudge=varargin{i+1};
            case 'CutOff'
                CutOff=varargin{i+1};
            case 'Points'
                Points=varargin{i+1};
            case 'ShowTailor'
                ShowTailor=varargin{i+1};
            case 'SelectByPortion'
                SelectByPortion=varargin{i+1};
            case 'Portion'
                Portion=varargin{i+1};
            case 'IfTailTailor'
                IfTailTailor=varargin{i+1};
            case 'IfFitExpTail'
                IfFitExpTail=varargin{i+1};
            case 'ExpTailPortion'
                ExpTailPortion=varargin{i+1};
            case 'Nsat'
                Nsat=varargin{i+1};
            case 'SM'
                SM=varargin{i+1};
            case 'IfBin'
                IfBin=varargin{i+1};
            case 'BinGridSize'
                BinGridSize=varargin{i+1};
            case 'IfSuperSampling'
                IfSuperSampling=varargin{i+1};
            case 'IfFourierFilter'
                IfFourierFilter=varargin{i+1};
            case 'CutOffFactor'
                CutOffFactor=varargin{i+1};
            case 'Nmethod'
                Nmethod=varargin{i+1};
            case 'pixellength'
                pixellength=varargin{i+1};
            case 'IfCleanImage'
                IfCleanImage=varargin{i+1};
            case 'OutlineExtrapolate'
                OutlineExtrapolate=varargin{i+1};
            case 'IfLookUpTable'
                IfLookUpTable=varargin{i+1};
            case 'ImagingTime'
                ImagingTime=varargin{i+1};
            case 'Zaveraging'
                Zaveraging=varargin{i+1};
            case 'omega'
                omega=varargin{i+1};
            case 'OutlineIntrapolate'
                OutlineIntrapolate=varargin{i+1};
            case 'IfExtrapolateAngle'
                IfExtrapolateAngle=varargin{i+1};
            case 'BoxShape'
                BoxShape=varargin{i+1};
        end
    end
end

%Defination of physical constant and experiment parameter
mli=9.988346*10^-27;  %kg
hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant


sigma0=0.215*10^-12/2; %in m^2


%Read the image
if ischar(Input)
    Img=fitsreadRL(Input);
    Nimg=AtomNumber( Img,pixellength^2,sigma0, Nsat);
    Nimg(isnan(Nimg))=0;
    Nimg(Nimg==inf)=0;
    Nimg(Nimg==-inf)=0;
else
    if size(Input,3)==1
        Nimg=Input;
    else
%         if ~IfLookUpTable
          Nimg=AtomNumber( Input,pixellength^2,sigma0, Nsat);
          IfLookUpTable=1;
%         else
%             
%             tic
%             Nimg=AtomNumberLUT( Input,pixellength^2,sigma0, Nsat, ImagingTime);
%             toc
%         end
    end
end
%Get the ROI1

if CropROI1
    questdlg('Now choose the ROI for get the n(z)')
    %msgbox('Now choose the ROI for get the n(z)');
    h=figure();
    imagesc(Nimg);
    caxis([0,15]);
    [~,Rect]=imcrop(h);
    close(h)
    ROI1=[round(Rect(1)),round(Rect(2)),round(Rect(1))+round(Rect(3)),round(Rect(2))+round(Rect(4))];
end
%Get the ROI2 for outline fitting
if CropROI2
    questdlg('Now choose the ROI for outline fitting')
    h=figure();
    imagesc(Nimg);
    caxis([0,15]);
    [~,Rect]=imcrop(h);
    close(h)
    ROI2=[round(Rect(1)),round(Rect(2)),round(Rect(1))+round(Rect(3)),round(Rect(2))+round(Rect(4))];
end
%ROI acquiring end

if IfLookUpTable && IfLookUpTable
    tic
    Nimg(ROI1(2):ROI1(4),ROI1(1):ROI1(3))=AtomNumberLUT( Input(ROI1(2):ROI1(4),ROI1(1):ROI1(3),:),pixellength^2,sigma0, Nsat, ImagingTime);
    toc
end


Nimg=Nimg-BGSubtraction;

%FourierFilter
if IfFourierFilter
    ImgFiltered=FourierFilter(Nimg(ROI1(2):ROI1(4),ROI1(1):ROI1(3)),CutOffFactor);
    Nimg(ROI1(2):ROI1(4),ROI1(1):ROI1(3))=ImgFiltered;
end



if IfCleanImage
    Nimg(ROI1(2):ROI1(4),ROI1(1):ROI1(3))=CleanImage(Nimg(ROI1(2):ROI1(4),ROI1(1):ROI1(3)));
end
Nimg(isnan(Nimg))=0;
Nimg(Nimg==inf)=0;
Nimg(Nimg==-inf)=0;


if IfSuperSampling
    Nimg=SuperSampling( Nimg,2,1)/4;
    pixellength=pixellength/2;
    ROI1=ROI1*2;
    ROI2=ROI2*2;
    zmin=zmin*2;
    zmax=zmax*2;
end

%Get the position of the tail, and tailor the tail to be a flat line
[n,z]=GenNvsZ( Nimg,ROI1,ROI2,pixellength,0,1 ,'ShowOutline',ShowOutline,'Nmethod',Nmethod,'OutlineExtrapolate',OutlineExtrapolate,'OutlineIntrapolate',OutlineIntrapolate,'Zaveraging',Zaveraging,'IfExtrapolateAngle',IfExtrapolateAngle,'BoxShape',BoxShape);


if CropTail
    h=figure();
    scatter(z,n);
    questdlg('Now give the range for tail fitting');
    [x,y]=getpts(h);
    close(h);
    zmin=min(x);
    zmax=max(x);
end

if IfTailTailor
    n=TailTailor(n,z,zmin,zmax,'ShowPlot',ShowTailor);
end

%Get the n(z)
n=fudge*n;
Ptf=TFfit(n,z,60);         
z0=Ptf(2); %center of the TF profile
Ztf=Ptf(3)*pixellength; %TF radius
ntf=TFfun( Ptf,z );
Z=(z-z0)*pixellength;

if IfBin
    Zmax=max(Z);
    Zgrid=linspace(0,Zmax,BinGridSize+1);
    [Z,n,~,~]=BinGrid( abs(Z),n,Zgrid,0);
    Z(isnan(n))=[];n(isnan(n))=[];
    Z=Z';n=n';
end

V=0.5*mli*omega^2*Z.^2;
Vtf=0.5*mli*omega^2*Ztf.^2;
output.Vtf=Vtf;

%excute the cutt off
Zcut=Ztf*CutOff;
n(abs(Z)>Zcut)=[];
V(abs(Z)>Zcut)=[];
z(abs(Z)>Zcut)=[];
ntf(abs(Z)>Zcut)=[];
Z(abs(Z)>Zcut)=[];

%sort Z,n with V
[Vsort,B]=sort(V);
nsort=n(B);
Zsort=Z(B);
zsort=z(B);

zcor.z=zsort;
zcor.z0=z0;

if IfHalf
    Vsort=Vsort(Zsort>0);
    nsort=nsort(Zsort>0);
    Zsort=Zsort(Zsort>0);
end

switch KappaMode
    case 0
        Kt=GetKappavsVPolyPoints(nsort,Vsort,Points,PolyOrder,'Smooth',ifsmooth);
    case 1
        Kt=GetKappavsVPolyVrange(nsort,Vsort,Vtf/VrangeFactor,PolyOrder,'Smooth',ifsmooth);
    case 2
        Kt=GetKappavsVPolyGlobal( nsort,Vsort,PolyOrder,'Smooth',ifsmooth );
    case 3
        Kt=GetKappavsVPolyWeighted(nsort,Vsort,Vtf/VrangeFactor,PolyOrder,'Smooth',ifsmooth);
    case 4
        Kt=GetKappavsVPolyExp( nsort,Vsort,0.66*Vtf,PolyOrder,'Smooth',ifsmooth);
    case 5
        Kt=GetKappavsVFiniteD( nsort,Vsort,SM,'Smooth',ifsmooth);
end

if IfFitExpTail
    [~,~,~,nf,Vf]=ExpFit_n( nsort,Vsort,ExpTailPortion );
else
    nf=nsort;
    Vf=Vsort;
end
[P,~,Pt]=GetPvsV( nf,Vf,'Smooth',ifsmooth );

%select the data points to show in Kt vs Pt data
%kmin=0.15;kmax=1.1;
nmax=max(nsort);
np=nsort/nmax;
if SelectByPortion
%     Ktsel=Kt(np>Portion);
%     Ptsel=Pt(np>Portion);
%     Zsel=Zsort(np>Portion);
%     Ktsel(abs(Zsel)<kmin*Ztf)=[];
%     Ptsel(abs(Zsel)<kmin*Ztf)=[];
%     Zsel(abs(Zsel)<kmin*Ztf)=[];
    mask1=np<Portion;
    mask2=abs(Zsort)>kmin*Ztf;
    mask=mask1 & mask2;
    Ktsel=Kt(mask);
    Ptsel=Pt(mask);
    Vsel=Vsort(mask);
else
%     Ktsel=Kt;Ptsel=Pt;Zsel=Zsort;
%     Ktsel(abs(Zsel)>kmax*Ztf)=[];
%     Ptsel(abs(Zsel)>kmax*Ztf)=[];
%     Zsel(abs(Zsel)>kmax*Ztf)=[];
%     Ktsel(abs(Zsel)<kmin*Ztf)=[];
%     Ptsel(abs(Zsel)<kmin*Ztf)=[];
%     Zsel(abs(Zsel)<kmin*Ztf)=[];
    mask1=abs(Zsort)>kmin*Ztf;
    mask2=abs(Zsort)<kmax*Ztf;
    mask=mask1 & mask2;
    Ktsel=Kt(mask);
    Ptsel=Pt(mask);
    Vsel=Vsort(mask);
end

%EFsort=real(hbar^2*(6*pi^2*nsort).^(2/3)/(2*mli));

EFsort=sign(nsort).*(hbar^2*(6*pi^2*abs(nsort)).^(2/3)/(2*mli));

%Plot nvsz with TF fitting
if FinalPlot
    
    h=figure();
    subplot(3,2,1);
    scatter(Z,n,'r.');
    hold on
    plot((z-z0)*pixellength,ntf,'b');
    hold off
    xlabel('z (pixel)');
    ylabel('n (m^{-3})');
    title('n vs z')
    %Plot n vs V
    subplot(3,2,2);
    scatter(Vsort/hh,EFsort/hh,'r.');
    xlabel('V (Hz)');
    ylabel('E_F (Hz)');
    title('EF vs V');
    %Plot Pt vs V
    subplot(3,2,3);
    scatter(Vsort/hh,Pt,'filled');
    ylim([0,5])
    xlabel('V (Hz)');
    ylabel('P/P_0 (m^{-3})');
    title('P/P_0 vs V');
    %Plot Kt vs V
    subplot(3,2,4);
    scatter(Vsort/hh,Kt,'filled');
    ylim([0,5])
    xlabel('V (Hz)');
    ylabel('\kappa/\kappa_0');
    title('\kappa/\kappa_0 vs V');
    %Plot Kt vs Pt
    subplot(3,2,5);
    scatter(Ptsel,Ktsel,'filled');
    ylim([0,5]);xlim([0,5]);
    xlabel('P/P_0');
    ylabel('\kappa/\kappa_0');
    title('EoS');
    %Plot P vs V
    subplot(3,2,6);
    scatter(Vsort/hh,P,'r.');
    ylabel('P');
    xlabel('V(Hz)');
    title('P vs V');
end
end

