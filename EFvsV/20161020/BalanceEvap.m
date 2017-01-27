GreenEvap=[1.60000000000000,1.40000000000000,1.20000000000000,1,0.800000000000000,0.600000000000000,0.400000000000000,0.200000000000000];
T=[1124,1074,906,772,655,533,340,150];
Terr=[80,85,52,92,78,23,21,15];
mu=[3951,3797,3581,3278,3324,2757,1922,968];
%%
T_mu=T./mu;
dT_mu=Terr./mu;

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
errorbar(GreenEvap,0.37*T_mu,0.37*dT_mu,'linestyle','none','marker','.','markersize',20)
xlim([0.0,1.8]);ylim([0,0.15])
xlabel('Green Evap (V)'); ylabel('0.37 \times k_B T / \mu')
set(axes1,'Ytick',[0,0.05,0.1,0.15,0.2])

%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
errorbar(GreenEvap,T,Terr,'linestyle','none','marker','.','markersize',20)
xlim([0.0,1.8]);
xlabel('Green Evap (V)'); ylabel('k_B T(Hz)')
%%
figure1 = figure;
axes1 = axes('Parent',figure1,'unit','inch','position',[1,1,4,3.2]);
plot(GreenEvap,mu,'linestyle','none','marker','.','markersize',20)
xlim([0.0,1.8]);
xlabel('Green Evap (V)'); ylabel('\mu(Hz)')