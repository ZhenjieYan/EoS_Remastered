GreenEvap=[1,1,0.800000000000000,0.600000000000000,1.20000000000000,1.40000000000000,1.60000000000000];
T=[496,292,318,287,469,513,596];
Terr=[31,83,44,55,72,60,76];
mu=[4709,5188,3943,3168,5345,5713,6300];
%%
T_mu=T./mu;
dT_mu=Terr./mu;

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
errorbar(GreenEvap,T_mu,dT_mu,'linestyle','none','marker','.','markersize',20)
xlim([0.4,1.8]);ylim([0,0.2])
xlabel('Green Evap (V)'); ylabel('k_B T / \mu')
set(axes1,'Ytick',[0,0.05,0.1,0.15,0.2])

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
errorbar(GreenEvap,T,Terr,'linestyle','none','marker','.','markersize',20)
xlim([0.4,1.8]);
xlabel('Green Evap (V)'); ylabel('k_B T(Hz)')
%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
plot(GreenEvap,mu,'linestyle','none','marker','.','markersize',20)
xlim([0.4,1.8]);
xlabel('Green Evap (V)'); ylabel('\mu(Hz)')