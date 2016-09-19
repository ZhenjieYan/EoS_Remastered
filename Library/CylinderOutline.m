function [x1,x2,X1,X2,Yt,p1,p2 ] = CylinderOutline( img,ROI,varargin )
% get the outline of the cylinder from top image
%img: the image need to fit
%ROI: [Xmin,Ymin,Xmax,Ymax] region where there is atom

Extrapolate=1;
Intrapolate=1;
IfExtrapolateAngle=0;

for i =1:length(varargin)
    if ischar(varargin{i})
    switch varargin{i}
        case 'Extrapolate'
            Extrapolate=varargin{i+1}; 
        case 'Intrapolate'
            Intrapolate=varargin{i+1}; 
        case 'IfExtrapolateAngle'
            IfExtrapolateAngle=varargin{i+1};
    end
    end
end


Rx = @(p,x) real(sqrt(p(2)^2- (x-p(1)).^2 ))*p(3);
% try
X=ROI(1):ROI(3);
Y=ROI(2)+1:ROI(4)-1;    
N=length(Y);
X1=zeros(1,N);X2=zeros(1,N); %X cordinator of the edge for each Y value
%imagesc(img);
%hold on

for i=1:N
    Nx=img(Y(i),X)+img(Y(i)+1,X)+img(Y(i)-1,X);
    P=OutlineFit(Nx,X,CMass1d(Nx,X),50);
    X1(i)=P(1)-abs(P(2));X2(i)=P(1)+abs(P(2));
    %plot(X1(i),Y(i),'r.','MarkerSize',20);
    %plot(X2(i),Y(i),'r.','MarkerSize',20);
end
Xc=(X1+X2)/2;
pc=polyfit(Y,Xc,1);


p1=polyfit(Y,X1,1);
p2=polyfit(Y,X2,1);
Yt=1:size(img,1);
xc=polyval(pc,Yt);

x1f=polyval(p1,Yt);
x2f=polyval(p2,Yt);
x1=x1f*0;
x2=x2f*0;

if Intrapolate
    x1((ROI(2)+1):(ROI(4)-1))=x1f((ROI(2)+1):(ROI(4)-1));
    x2((ROI(2)+1):(ROI(4)-1))=x2f((ROI(2)+1):(ROI(4)-1));
else
    x1((ROI(2)+1):(ROI(4)-1))=X1;
    x2((ROI(2)+1):(ROI(4)-1))=X2;
end

if (~Extrapolate)
    Rup=(x2f(Yt==(ROI(4)-1))-x1f(Yt==(ROI(4)-1)))/2;
    Rdown=(x2f(Yt==(ROI(2)+1))-x1f(Yt==(ROI(2)+1)))/2;
    if ~IfExtrapolateAngle
        x1(Yt<=ROI(2))=xc(Yt==(ROI(2)+1))-Rdown;
        x2(Yt<=ROI(2))=xc(Yt==(ROI(2)+1))+Rdown;
        x1(Yt>=ROI(4))=xc(Yt==(ROI(4)-1))-Rup;
        x2(Yt>=ROI(4))=xc(Yt==(ROI(4)-1))+Rup;
    else
        x1(Yt<=ROI(2))=xc(Yt<=ROI(2))-Rdown;
        x2(Yt<=ROI(2))=xc(Yt<=ROI(2))+Rdown;
        x1(Yt>=ROI(4))=xc(Yt>=ROI(4))-Rup;
        x2(Yt>=ROI(4))=xc(Yt>=ROI(4))+Rup;
    end
%     x1(Yt<=ROI(2))=x1f(Yt==(ROI(2)+1));
%     x1(Yt>=ROI(4))=x1f(Yt==(ROI(4)-1));
%     x2(Yt<=ROI(2))=x2f(Yt==(ROI(2)+1));
%     x2(Yt>=ROI(4))=x2f(Yt==(ROI(4)-1));
else
    x1(Yt<=ROI(2))=x1f(Yt<=ROI(2));
    x1(Yt>=ROI(4))=x1f(Yt>=ROI(4));
    x2(Yt<=ROI(2))=x2f(Yt<=ROI(2));
    x2(Yt>=ROI(4))=x2f(Yt>=ROI(4));
end

end

