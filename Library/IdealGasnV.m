function n = IdealGasnV( P,V )
% P(1): Chemical potential of the gas;
% P(2): k_B*Temperature of the gas
% Physical Constant:
mli=9.988346*10^-27;  %kg
hbar=1.0545718*10^(-34); %SI
hh=2*pi*hbar;%SI Planck constant

mu_local=P(1)-V;
beta_mu=mu_local/P(2);
Z_vec = exp(beta_mu);

TTilde = (4*pi)./(6*pi^2*(-PolyLogFrac(3/2,-Z_vec))).^(2/3);

EF=P(2)./TTilde;

n=(2*mli*EF*hh/hbar^2).^(3/2)/(6*pi^2)/1e18;
% mask1=(n==-inf);
% mask2=(n==inf);
% mask3=isnan(n);
% 
% mask=mask1 & mask2 & mask3;
% 
% if sum(mask)>0
%     disp('????')
% end
% 
% n(mask)=0;

end

