function P = IdealGasFitwoFudge( V,n,P0 )
% P(1): Chemical potential of the gas;
% P(2): k_B*Temperature of the gas

FitFun=@(beta,x) IdealGasnV( [beta(1),beta(2)],x );

P = nlinfit(V,n,FitFun,P0);


end

