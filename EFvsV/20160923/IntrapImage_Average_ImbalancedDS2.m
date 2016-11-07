%%

filefolder_Unitary='/Users/Zhenjie/Data/2016-02-23/';
load('/Users/Zhenjie/Data/Processed/2016-09-23/2016-09-23Set2Bin.mat')
addpath('../../Library');
addpath('/Users/Zhenjie/Github/BoxPaper/NewHybrid/Code/ColorMap');
%% State1
NimgAverage=0;
for i=1:length(imglistS1Final)
    RawImg=imglistS1Final{i};
    sigma0=0.215*10^-12/2; %in m^2
    Nsat=770;
    pixellength=0.7*3*10^-6; %in m
    Nimg=AtomNumberLUT( RawImg,pixellength^2,sigma0, Nsat,10)*1.7087;
    NimgAverage=NimgAverage+Nimg;
end
BGimg=0;
NimgAverage=NimgAverage/length(imglistS1Final);
NimgAverage=NimgAverage-BGimg;
%NimgAverage=CleanImage(NimgAverage);
%% line of sight: get outline
ROI2=[320,380,456,500];
[x1,x2,X1,X2,Yt,p1,p2 ] = CylinderOutline( NimgAverage,ROI2 );
imagesc(NimgAverage);
hold on
plot(x1,Yt,'r.','markersize',5);
plot(x2,Yt,'r.','markersize',5);
hold off
caxis([0,120]);
Xc=(x1+x2)/2;
R=(x2-x1)/2;
%% line of sight: reconstruction state1
Ncor=NimgAverage*0;
for i=1:size(NimgAverage,1)
    for j=(ceil(x1(i))+0):(floor(x2(i))-0)
        Ri=R(i);
        hj=j-Xc(i);
        Zij=sqrt(Ri^2-hj^2);
        Ncor(i,j)=NimgAverage(i,j)/Zij;
    end
    disp(i)
end

imagesc(Ncor)
caxis([-0.4,3])

%% plot the un-correct version state1

ROI=[300,265,456,600];
z0=430.84;
x0=385.5;
Z=((ROI(2):ROI(4))-z0)*2.1;
X=((ROI(1):ROI(3))-x0)*2.1;


pic=NimgAverage(ROI(2):ROI(4),ROI(1):ROI(3));
imagesc(pic)

pic=pic';

% Create figure
figure1 = figure;
%colormap(colormap_boxpaper(100));

% Create axes
axes1 = axes('Parent',figure1)%,'unit','inch','position',[1,1,1.4,1.4]);

% Create image
image(Z,X,pic,'Parent',axes1,'CDataMapping','scaled');
hold on

box(axes1,'on');

set(axes1,'CLim',[0 1],'DataAspectRatio',[1 1 1],'Layer','top','Ydir','normal');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
caxis([-5,120]);
% ylim([-100,100]);
% xlim([-200,200])
hold off
%% Shear transform state1
tform = maketform('affine',[1 0 0; 0.12 1 0; 0 0 1]);
ImgTran= imtransform(Ncor,tform);
imagesc(ImgTran);


%% plot the line of sight corrected version state1

ROI=[385,292,500,574];
z0=430.84;
x0=438;
Z=((ROI(2):ROI(4))-z0)*2.1;
X=((ROI(1):ROI(3))-x0)*2.1;


pic=ImgTran(ROI(2):ROI(4),ROI(1):ROI(3));
imagesc(pic)

pic=pic';
colormap_Unitary=cbrewer('seq', 'Blues', 200, 'cubic');
% Create figure
figure1 = figure;
colormap(colormap_Unitary);

% Create axes
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,1.47]);

% Create image
image(Z,X,pic,'Parent',axes1,'CDataMapping','scaled');
hold on
% line([20,20+100/1.4],[30,30],'linewidth',3,'color',[1,1,1])
% line([20,20],[22,38],'linewidth',1,'color',[1,1,1]);
% line([20+100/1.4,20+100/1.4],[22,38],'linewidth',1,'color',[1,1,1])
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 281.5]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.5 426.5]);
box(axes1,'on');
% line([Zcr1,Zcr1],[-1000,1000],'color',[201,67,52]/255,'linewidth',1);
% line([Zcr2,Zcr2],[-1000,1000],'color',[201,67,52]/255,'linewidth',1)
%xlim([1,size(pic,2)]);ylim([1,size(pic,1)])
% Set the remaining axes properties
set(axes1,'CLim',[0 1],'DataAspectRatio',[1 1 1],'Layer','top','Ydir','normal');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
caxis([-0.2,5]);
ylim([0,50]);
xlim([-250,250])
hold off
%%
save('DS2Maj.mat','Z','X','pic');

%% State 2

NimgAverage=0;
for i=1:length(imglistS2Final)
    RawImg=imglistS2Final{i};
    sigma0=0.215*10^-12/2; %in m^2
    Nsat=770;
    pixellength=0.7*3*10^-6; %in m
    Nimg=AtomNumberLUT( RawImg,pixellength^2,sigma0, Nsat,10)*1.7087;
    NimgAverage=NimgAverage+Nimg;
end
BGimg=0;
NimgAverage=NimgAverage/length(imglistS2Final);
NimgAverage=NimgAverage-BGimg;
%NimgAverage=CleanImage(NimgAverage);
%% line of sight: get outline
ROI2=[270,428,456,455];
[x1,x2,X1,X2,Yt,p1,p2 ] = CylinderOutline( NimgAverage,ROI2 );
imagesc(NimgAverage);
hold on
plot(x1,Yt,'r.','markersize',5);
plot(x2,Yt,'r.','markersize',5);
hold off
caxis([0,120]);
Xc=(x1+x2)/2;
R=(x2-x1)/2;
%% line of sight: reconstruction state1
Ncor=NimgAverage*0;
for i=1:size(NimgAverage,1)
    for j=(ceil(x1(i))+0):(floor(x2(i))-0)
        Ri=R(i);
        hj=j-Xc(i);
        Zij=sqrt(Ri^2-hj^2);
        Ncor(i,j)=NimgAverage(i,j)/Zij;
    end
    disp(i)
end

imagesc(Ncor)
caxis([-0.4,3])

%% plot the un-correct version state1

ROI=[300,265,456,600];
z0=430.84;
x0=385.5;
Z=((ROI(2):ROI(4))-z0)*2.1;
X=((ROI(1):ROI(3))-x0)*2.1;


pic=NimgAverage(ROI(2):ROI(4),ROI(1):ROI(3));
imagesc(pic)

pic=pic';

% Create figure
figure1 = figure;
%colormap(colormap_boxpaper(100));

% Create axes
axes1 = axes('Parent',figure1)%,'unit','inch','position',[1,1,1.4,1.4]);

% Create image
image(Z,X,pic,'Parent',axes1,'CDataMapping','scaled');
hold on

box(axes1,'on');

set(axes1,'CLim',[0 1],'DataAspectRatio',[1 1 1],'Layer','top','Ydir','normal');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
caxis([-5,120]);
% ylim([-100,100]);
% xlim([-200,200])
hold off
%% Shear transform state1
tform = maketform('affine',[1 0 0; 0.12 1 0; 0 0 1]);
ImgTran= imtransform(Ncor,tform);
imagesc(ImgTran);


%% plot the line of sight corrected version state1

ROI=[385,292,500,574];
z0=430.84;
x0=438;
Z=((ROI(2):ROI(4))-z0)*2.1;
X=((ROI(1):ROI(3))-x0)*2.1;


pic=ImgTran(ROI(2):ROI(4),ROI(1):ROI(3));
imagesc(pic)

pic=pic';
colormap_Unitary=cbrewer('seq', 'Blues', 200, 'cubic');

temp1=colormap_Unitary(:,1);
temp2=colormap_Unitary(:,2);
temp3=colormap_Unitary(:,3);

colormap_Unitary(:,3)=temp1;
colormap_Unitary(:,1)=temp3;
colormap_Unitary(:,2)=temp1;

% Create figure
figure1 = figure;
colormap(colormap_Unitary);

% Create axes
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,1.47]);

% Create image
image(Z,X,pic,'Parent',axes1,'CDataMapping','scaled');
hold on
% line([20,20+100/1.4],[30,30],'linewidth',3,'color',[1,1,1])
% line([20,20],[22,38],'linewidth',1,'color',[1,1,1]);
% line([20+100/1.4,20+100/1.4],[22,38],'linewidth',1,'color',[1,1,1])
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[0.5 281.5]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[0.5 426.5]);
box(axes1,'on');
% line([Zcr1,Zcr1],[-1000,1000],'color',[201,67,52]/255,'linewidth',1);
% line([Zcr2,Zcr2],[-1000,1000],'color',[201,67,52]/255,'linewidth',1)
%xlim([1,size(pic,2)]);ylim([1,size(pic,1)])
% Set the remaining axes properties
set(axes1,'CLim',[0 1],'DataAspectRatio',[1 1 1],'Layer','top','Ydir','normal');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
caxis([0,5]);
ylim([-50,0]);
xlim([-250,250])
hold off
%%
save('DS2Min.mat','Z','X','pic')