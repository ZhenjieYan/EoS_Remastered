function imgout = ImgRotate( img,angle )
% Do the roatetion for the raw image
% img: raw image;
% angle: Rotation angle in degrees;

m=size(img,3);
temp=imrotate(img(:,:,1),angle);
imgout=zeros(size(temp,1),size(temp,2),m);

for i=1:m
    imgout(:,:,i)=imrotate(img(:,:,i),angle,'bilinear');
end


end

