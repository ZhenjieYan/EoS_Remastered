function P=OutlineFitSquare(Nx,x,x0,R0)
%fitfun = @(p,x) real(sqrt(p(2)^2- (x-p(1)).^2 ))*p(3);
squarefun=@(P,y) 0.5*(erf((y-P(1)+P(2))/P(4))-erf((y-P(1)-P(2))/P(4)))*P(3);
P=nlinfit(x,Nx,squarefun,[x0,R0,max(Nx),2]);
end