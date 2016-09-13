function Nimg = AtomNumberLUT(img,pixelsize,sigma, Nsat,time  )
%ATOMNUMBERLUT Summary of this function goes here
%   Detailed explanation goes here
DictionaryFolder='/Users/Zhenjie/Data/LookUpTable/';

load([DictionaryFolder,'DictionaryT',num2str(time,'%02.0f'),'.mat']);

WA=(img(:,:,1)-img(:,:,3))/Nsat;
WOA=(img(:,:,2)-img(:,:,3))/Nsat;

ODtemp=griddata(Iini,Ifin,OD,WOA,WA);

Nimg=ODtemp/sigma*pixelsize;

mask1=isnan(Nimg);
mask2=WOA<max(Iini);
mask3=(WA>0) & (WOA>0);
mask4=(WOA<=WA);

mask=mask1&mask2&mask3&mask4;
Nimg(mask)=-griddata(Iini,Ifin,OD,WA(mask),WOA(mask));


mask1=isnan(Nimg);
mask2=(WA>0) & (WOA>0);
mask3=WOA<max(Iini);
mask4=((WA>0.95*WOA)&(WA<1.05*WOA)) | ((WOA<0.05)&(WA<0.05));
mask=mask1&mask2&mask3&mask4;

OD_con=log(WOA./WA);
Nimg(mask)=OD_con(mask);

end

