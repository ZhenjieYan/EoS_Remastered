addpath('C:\Users\BEC1\Documents\GitHub\EoS_Remastered\Library')
%%
mli=9.988346*10^-27; %kg
hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant
omega=23.9*2*pi; %in rad/s
pixellength=0.7*10^-6*3; %in m
sigma0=0.215*10^-12/2; %in m^2
%%
imagelist={};
for i=15:22
    load(['/Users/Zhenjie/Data/Processed/2016-10-20/2016-10-20Set',num2str(i),'Bin.mat']);
    imglisttemp=imglistS1Final;
    imagelist=[imagelist,imglisttemp];
    disp(i)
end
%%
Vtf=zeros(length(imagelist),1);
R=Vtf;
outputlist=cell(length(imagelist),1);
for i=1:length(imagelist)
    [Pt,Kt,nsort,Vsort,Zsort,Ptsel,Ktsel,EFS1,P,zcor,Vsel,output]=EOS_Online(imagelist{i},'ROI1',[220,265,426,630],...
    'ROI2',[270,430,456,455],'ShowOutline',0,'TailRange',[325,575],'KappaMode',2,'PolyOrder',10,'VrangeFactor',5,'IfHalf',0,'kmax',0.9,'kmin',0.15,...
    'Fudge',1.25,'BGSubtraction',0,'IfFitExpTail',0,'Nsat',770,'ShowPlot',0,'CutOff',inf,'IfHalf',0,'pixellength',0.7e-6*3,'SM',3,'IfBin',0,'BinGridSize',150,...
    'IfCleanImage',1,'OutlineExtrapolate',1,'IfLookUpTable',1,'Zaveraging',0,'BoxShape','Square');
    Vtf(i)=output.Vtf;
    outputlist{i}=output;
    disp(i)
end
Vtf=Vtf/hh;
%%
outlinelist=cell(length(imagelist),1);
D=zeros(length(imagelist),1);
YI=448;
for i=1:length(outputlist)
    outlinelist{i}=outputlist{i}.outline;
    D(i)=outputlist{i}.outline.x2(YI)-outputlist{i}.outline.x1(YI);
end