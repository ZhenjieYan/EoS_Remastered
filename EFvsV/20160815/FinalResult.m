DataList=dir('/Users/Zhenjie/Github/EoS_Remastered/EFvsV/20160815/Processed');
addpath('../../Library');

NameList={};

KappaTildeS2List=[];
nS2List=[];
nS1List=[];
TTildeS1List=[];
SigmaList=[];

s=0;

for i=1:length(DataList)
    Datatemp=DataList(i);
    Nametemp=Datatemp.name;
    if length(Nametemp)>3
        if strcmp(Nametemp(end-2:end),'mat')
            NameList=[NameList;Nametemp];
        end
    end
end

Th=1;
for i=1:length(NameList)
    Datatemp=load(['Processed/',NameList{i}]);
    KappaTildeS2List=[KappaTildeS2List,Datatemp.KappaTildeS2BinV(Th:end)];
    nS2List=[nS2List,Datatemp.nS2BinV(Th:end)];
    nS1List=[nS1List,Datatemp.nS1BinV(Th:end)];
    TTildeS1List=[TTildeS1List,Datatemp.T_trap./Datatemp.EFS1BinV(Th:end)];
end

%%
mask1=nS1List>0.24;
mask2=nS1List<0.005;
mask=mask1 | mask2;

KappaTildeS2List(mask)=[];
nS2List(mask)=[];
nS1List(mask)=[];
TTildeS1List(mask)=[];

SigmaList=(nS1List-nS2List)./(nS1List+nS2List);
%%
scatter3(SigmaList,TTildeS1List,KappaTildeS2List);
xlim([0,1]);ylim([0,0.3]);
xlabel('\sigma');ylabel('T/T_{F,\uparrow}')
%%
mask1=SigmaList<1;
mask2=SigmaList>0;
mask3=TTildeS1List<0.3;
mask4=TTildeS1List>0;
mask=mask1&mask2&mask3&mask4;

SigmaShow=SigmaList(mask);
TTildeS1Show=TTildeS1List(mask);
KappaTildeS2Show=KappaTildeS2List(mask);


scatter3(SigmaShow,TTildeS1Show,KappaTildeS2Show);
xlim([0,1]);ylim([0,0.3]);
xlabel('\sigma');ylabel('T/T_{F,\uparrow}')
%%
[KappaTildeS2grid,SigmaFit,TTildeS1Fit] = gridfit(SigmaShow,TTildeS1Show,KappaTildeS2Show,100,100);
surf(SigmaFit,TTildeS1Fit,KappaTildeS2grid);
xlim([0,1]);ylim([0,0.3]);
xlabel('\sigma');ylabel('T/T_{F,\uparrow}')

%%
h=figure();
axes1=axes('parent',h);
imagesc(SigmaFit(1,:),TTildeS1Fit(:,1),KappaTildeS2grid);
caxis([0,4])
k=1.08;
hold on
plot(SigmaShow,TTildeS1Show,'k.')
line([0,0.39],[0.167,0],'linewidth',2);
line([0,0.2]*k,[0,0.07]*k,'linewidth',2);
text(0.01,0.073,'Superfluid','color','b','FontSize',14);
text(0.13,0.03,'Unstable','color',[62,166,60]/255,'FontSize',14);
text(0.4,0.08,'Normalfluid','color',[201,67,52]/255,'FontSize',14);
hold off
set(axes1,'YDir','normal');
colorbar()
xlim([0,1]);ylim([0,0.3]);
xlabel('(n_1-n_2)/(n_1+n_2)');
ylabel('T/T_{F,1}');

title('Compressibility of the minority')

%%
Th=1
DataSet='Result0815Set19.mat';
Datafocus=load(['Processed/',DataSet]);
KappaTildeS2focus=Datafocus.KappaTildeS2BinV(Th:end);
nS2focus=Datafocus.nS2BinV(Th:end);
nS1focus=Datafocus.nS1BinV(Th:end);
TTildeS1focus=Datafocus.T_trap./Datafocus.EFS1BinV(Th:end);
Sigmafocus=(nS1focus-nS2focus)./(nS1focus+nS2focus);

T_focus=Datafocus.T_trap
mu_focus=Datafocus.mu_trap_S1
%%
h=figure();
axes1=axes('parent',h);
imagesc(SigmaFit(1,:),TTildeS1Fit(:,1),KappaTildeS2grid);
caxis([0,4])
k=1.08;
hold on
plot(SigmaShow,TTildeS1Show,'k.');
line([0,0.39],[0.167,0],'linewidth',2);
line([0,0.2]*k,[0,0.07]*k,'linewidth',2);
text(0.01,0.073,'Superfluid','color','b','FontSize',14);
text(0.13,0.03,'Unstable','color',[62,166,60]/255,'FontSize',14);
text(0.4,0.08,'Normalfluid','color',[201,67,52]/255,'FontSize',14);
plot(Sigmafocus,TTildeS1focus,'ro','linewidth',2);
hold off
set(axes1,'YDir','normal');
colorbar()
xlim([0,1]);ylim([0,0.3]);
xlabel('(n_1-n_2)/(n_1+n_2)');
ylabel('T/T_{F,1}');

title('Compressibility of the minority')