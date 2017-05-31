load('DS2NonSliced.mat');
ZBinZ_NonSliced=ZBinZ;
nS1BinZ_NonSliced=nS1BinZ;
nS2BinZ_NonSliced=nS2BinZ;
%%
load('DS2Sliced.mat');
ZBinZ_Sliced=ZBinZ;
nS1BinZ_Sliced=nS1BinZ;
nS2BinZ_Sliced=nS2BinZ;
%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
plot(ZBinZ_NonSliced,nS1BinZ_NonSliced,'DisplayName','State1, Before Depumping','linewidth',1.5);
hold on
plot(ZBinZ_Sliced,nS1BinZ_Sliced,'DisplayName','State1, After Depumping','linewidth',1.5);
hold off
grid on
xlabel('Z(um)');ylabel('n(z) (um^{-3})');
legend show
ylim([-0.05,0.27])
%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
plot(ZBinZ_NonSliced,nS2BinZ_NonSliced,'DisplayName','State2, Before Depumping','linewidth',1.5);
hold on
plot(ZBinZ_Sliced,nS2BinZ_Sliced,'DisplayName','State2, After Depumping','linewidth',1.5);
hold off
grid on
xlabel('Z(um)');ylabel('n(z) (um^{-3})');
legend show
ylim([-0.05,0.27])
%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
ZBinZ_com=ZBinZ_NonSliced;
nS1BinZ_Sliced_com=interp1(ZBinZ_Sliced,nS1BinZ_Sliced,ZBinZ_com);
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);
plot(ZBinZ_com,nS1BinZ_Sliced_com./nS1BinZ_NonSliced,'DisplayName','State1, Before Depumping','linewidth',1.5);
ylim([0.6,0.8]);xlim([-100,100]);
xlabel('Z(um)');ylabel('n_{after}/n_{before}');

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
ZBinZ_com=ZBinZ_NonSliced;
nS2BinZ_Sliced_com=interp1(ZBinZ_Sliced,nS2BinZ_Sliced,ZBinZ_com);
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,3,1.8]);
plot(ZBinZ_com,nS2BinZ_Sliced_com./nS2BinZ_NonSliced,'DisplayName','State1, Before Depumping','linewidth',1.5);
ylim([0.8,1]);xlim([-100,100]);
xlabel('Z(um)');ylabel('n_{after}/n_{before}');