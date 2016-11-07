%%
addpath('../../Library');
addpath('/Users/Zhenjie/Github/BoxPaper/NewHybrid/Code/ColorMap');
Power=1/2;
%% plot the line of sight corrected version state2
load('DS2Min.mat')
Imgstrip=pic(abs(X)<40,:);
Nz=sum(Imgstrip)/size(Imgstrip,1);
Nbg=mean(Nz(abs(Z)>150));
pic=pic-Nbg;

%%
picpPower=abs(pic.^Power).*sign(pic);
ImgPowerstrip=picpPower(abs(X)<40,:);
NzPower=sum(ImgPowerstrip)/size(ImgPowerstrip,1);
plot(Z,NzPower)
NbgPower=mean(NzPower(abs(Z)>150));
picpPower=picpPower-NbgPower;
%%
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
image(Z,X,picpPower,'Parent',axes1,'CDataMapping','scaled');
hold on
box(axes1,'on');
set(axes1,'CLim',[0 1],'DataAspectRatio',[1 1 1],'Layer','top','Ydir','normal');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
caxis([0,6^Power]);
ylim([-50,0]);
xlim([-250,250])
hold off

%% plot the line of sight corrected version state1
load('DS2Maj.mat')
Imgstrip=pic(abs(X)<40,:);
Nz=sum(Imgstrip)/size(Imgstrip,1);
Nbg=mean(Nz(abs(Z)>200));
pic=pic-Nbg;
%%
picpPower=abs(pic.^Power).*sign(pic);
ImgPowerstrip=picpPower(abs(X)<40,:);
NzPower=sum(ImgPowerstrip)/size(ImgPowerstrip,1);
plot(Z,NzPower)
NbgPower=mean(NzPower(abs(Z)>150));
picpPower=picpPower-NbgPower;


%%
colormap_Unitary=cbrewer('seq', 'Blues', 200, 'cubic');

% temp1=colormap_Unitary(:,1);
% temp2=colormap_Unitary(:,2);
% temp3=colormap_Unitary(:,3);
% 
% colormap_Unitary(:,3)=temp1;
% colormap_Unitary(:,1)=temp3;
% colormap_Unitary(:,2)=temp1;

% Create figure
figure1 = figure;
colormap(colormap_Unitary);

% Create axes
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,1.47,1.47]);

% Create image
image(Z,X,picpPower,'Parent',axes1,'CDataMapping','scaled');
hold on
box(axes1,'on');
set(axes1,'CLim',[0 1],'DataAspectRatio',[1 1 1],'Layer','top','Ydir','normal');
set(axes1,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
caxis([0,5^Power]);
ylim([-50,0]);
xlim([-250,250])
hold off
