% Calculate PSNR
% off: offset (default: 5)
% maxval: maximum value (default: 255)

function psnr = psnrfun(image1, image2, off, maxval);

if nargin<4
    maxval = 255;
end
if nargin<3
    off = 5;
end
  
image1 = double(image1);
s = size(image1);
image1 = image1(off(1)+1:(s(1)-off(1)), off(2)+1:(s(2)-off(2)));


image2 = double(image2);
image2 = image2(off(1)+1:(s(1)-off(1)), off(2)+1:(s(2)-off(2)));

x = image1 - image2;
x = x.^2;
s1=size(image1);
y = sum(sum(x))/(s1(1)*s1(2));
y = y.^0.5;
psnr = 20*log10(maxval/y);                                           
