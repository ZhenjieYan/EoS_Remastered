function ImageOut = CleanImage( Img )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ImageOut=Img;
[Yrange,Xrange]=size(Img);
[X,Y]=meshgrid(1:Xrange,1:Yrange);

ImgColumn=Img(:);
XColumn=X(:);
YColumn=Y(:);

Xabnormal_list=[];
Yabnormal_list=[];

Xabnormal_list=[Xabnormal_list;XColumn(ImgColumn==inf)];
Xabnormal_list=[Xabnormal_list;XColumn(ImgColumn==-inf)];
Xabnormal_list=[Xabnormal_list;XColumn(isnan(ImgColumn))];

Yabnormal_list=[Yabnormal_list;YColumn(ImgColumn==inf)];
Yabnormal_list=[Yabnormal_list;YColumn(ImgColumn==-inf)];
Yabnormal_list=[Yabnormal_list;YColumn(isnan(ImgColumn))];

N=length(Xabnormal_list);
for i=1:N
    Xi=Xabnormal_list(i);
    Yi=Yabnormal_list(i);
    
    mask1=abs(XColumn-Xi)<2;
    mask2=abs(YColumn-Yi)<2;
    mask=mask1 & mask2;
    
    Neighbor=ImgColumn(mask);
    
    Neighbor(isnan(Neighbor))=[];
    Neighbor(Neighbor==inf)=[];
    Neighbor(Neighbor==-inf)=[];
    r=mean(Neighbor);
    if isnan(r)
        r=0;
    end
    ImageOut(Yi,Xi)=r;
end


end

