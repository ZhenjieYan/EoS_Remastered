function imgout = ImgBin( img,Bin )
% Do the Binning for the raw image
% img: raw image;
% Bin: Bin factor;

m=size(img,3);
temp=imresize(img(:,:,1),1/Bin,'Box');
imgout=zeros(size(temp,1),size(temp,2),m);
for i=1:m
    imgout(:,:,i)=imresize(img(:,:,i),1/Bin,'Box');
end

end

