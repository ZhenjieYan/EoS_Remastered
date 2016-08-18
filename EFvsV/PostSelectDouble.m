%%
filelist={'08-15-2016_16_32_12_TopA';'08-15-2016_16_32_12_TopB';'08-15-2016_16_31_19_TopA';'08-15-2016_16_31_19_TopB';'08-15-2016_16_29_34_TopA';'08-15-2016_16_29_34_TopB';'08-15-2016_16_28_42_TopA';'08-15-2016_16_28_42_TopB';'08-15-2016_16_27_49_TopA';'08-15-2016_16_27_49_TopB';'08-15-2016_16_26_05_TopA';'08-15-2016_16_26_05_TopB';'08-15-2016_16_25_12_TopA';'08-15-2016_16_25_12_TopB';'08-15-2016_16_24_20_TopA';'08-15-2016_16_24_20_TopB';'08-15-2016_16_23_27_TopA';'08-15-2016_16_23_27_TopB'};
folder='/Users/Zhenjie/Data/2016-08-15/';
%%
%filelist=filelist(1:2);
N=length(filelist);
filelistS1=filelist(1:2:N-1);
filelistS2=filelist(2:2:N);
%%
%Define the physical constant
mli=9.988346*10^-27;  %kg
hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant
omega=23.9*2*pi; %in rad/s
pixellength=0.7*10^-6; %in m
sigma0=0.215*10^-12/2; %in m^2
%load all the functions
addpath('../Library');
Nsat=120;

%% Get the total atom number for S2 img;

ROI=[850,550,1300,1300];
Z=ROI(2):ROI(4);
Z=Z';
zmin=650;zmax=1200;
indexS2=1:length(filelistS2);
imglistS2={};
numlistS2=[];
for i=1:length(filelistS2)
    tempimg=fitsread([folder,filelistS2{i},'.fits']);
    imglistS2=[imglistS2,tempimg];
    Nimg=AtomNumber(tempimg,pixellength^2,sigma0, Nsat);
    Nimg=Nimg(ROI(2):ROI(4),ROI(1):ROI(3));
    Nimg=CleanImage(Nimg);
    Nz=sum(Nimg,2);
    Nz=TailTailor(Nz,Z,zmin,zmax);
    numlistS2=[numlistS2;sum(Nz)];
    disp(i);
end
scatter(indexS2,numlistS2)
ylim([0,1.1*max(numlistS2)]);
xlabel('Img index');ylabel('Total Atom number');

%%
tolerance=0.05;
minscan=min(numlistS2);maxscan=max(numlistS2);
MeanNumlist=linspace(minscan,maxscan,40);
Npick=0*MeanNumlist;

for i=1:length(MeanNumlist)
    mask1=numlistS2<=(MeanNumlist(i)*(1+tolerance));
    mask2=numlistS2>=(MeanNumlist(i)*(1-tolerance));
    mask=mask1 & mask2;
    Npick(i)=sum(mask);
end
scatter(MeanNumlist,Npick);

[~,B]=max(Npick);
MeanNum=MeanNumlist(B);
mask1=numlistS2<=(MeanNum*(1+tolerance));
mask2=numlistS2>=(MeanNum*(1-tolerance));
mask=mask1 & mask2;

imglistS2pick1=imglistS2(mask);
filelistS2pick1=filelistS2(mask);
filelistS1pick1=filelistS1(mask);

%% Check the atom number for state 1

ROI=[850,350,1300,1300];
Z=ROI(2):ROI(4);
Z=Z';
zmin=500;zmax=1200;
indexS1=1:length(filelistS1pick1);
imglistS1={};
numlistS1=[];
for i=1:length(filelistS1pick1)
    tempimg=fitsread([folder,filelistS1pick1{i},'.fits']);
    imglistS1=[imglistS1,tempimg];
    Nimg=AtomNumber(tempimg,pixellength^2,sigma0, Nsat);
    Nimg=Nimg(ROI(2):ROI(4),ROI(1):ROI(3));
    Nimg=CleanImage(Nimg);
    Nz=sum(Nimg,2);
    Nz=TailTailor(Nz,Z,zmin,zmax);
    numlistS1=[numlistS1;sum(Nz)];
    disp(i);
end
scatter(indexS1,numlistS1)
ylim([0,1.1*max(numlistS1)]);
xlabel('Img index');ylabel('Total Atom number');
%%
tolerance=0.05;
minscan=min(numlistS1);maxscan=max(numlistS1);
MeanNumlist=linspace(minscan,maxscan,40);
Npick=0*MeanNumlist;

for i=1:length(MeanNumlist)
    mask1=numlistS1<=(MeanNumlist(i)*(1+tolerance));
    mask2=numlistS1>=(MeanNumlist(i)*(1-tolerance));
    mask=mask1 & mask2;
    Npick(i)=sum(mask);
end
scatter(MeanNumlist,Npick);

[~,B]=max(Npick);
MeanNum=MeanNumlist(B);
mask1=numlistS1<=(MeanNum*(1+tolerance));
mask2=numlistS1>=(MeanNum*(1-tolerance));
mask=mask1 & mask2;

filelistS1pick2=filelistS1pick1(mask);
filelistS2pick2=filelistS2pick1(mask);
%%
ROI=[850,350,1300,1300];
Z=ROI(2):ROI(4);
Z=Z';
zmin=500;zmax=1200;
Index=1:length(filelistS2pick2);
imglistS1Final={};
numlistS1Final=[];
imglistS2Final={};
numlistS2Final=[];

for i=1:length(filelistS1pick2)
    tempimgS1=fitsread([folder,filelistS1pick2{i},'.fits']);
    imglistS1Final=[imglistS1Final,tempimgS1];
    Nimg=AtomNumber(tempimgS1,pixellength^2,sigma0, Nsat);
    Nimg=Nimg(ROI(2):ROI(4),ROI(1):ROI(3));
    Nimg=CleanImage(Nimg);
    Nz=sum(Nimg,2);
    Nz=TailTailor(Nz,Z,zmin,zmax);
    numlistS1Final=[numlistS1Final;sum(Nz)];
    
    tempimgS2=fitsread([folder,filelistS2pick2{i},'.fits']);
    imglistS2Final=[imglistS2Final,tempimgS2];
    Nimg=AtomNumber(tempimgS2,pixellength^2,sigma0, Nsat);
    Nimg=Nimg(ROI(2):ROI(4),ROI(1):ROI(3));
    Nimg=CleanImage(Nimg);
    Nz=sum(Nimg,2);
    Nz=TailTailor(Nz,Z,zmin,zmax);
    numlistS2Final=[numlistS2Final;sum(Nz)];
    
    disp(i);
end
%%
scatter(Index,numlistS1Final);
hold on
scatter(Index,numlistS2Final);
hold off

filelistS1Final=filelistS1pick2;
filelistS2Final=filelistS2pick2;
%%
save('2016-08-15set2.mat','imglistS1Final','imglistS2Final','filelistS1Final','filelistS2Final')
