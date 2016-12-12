%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function im_new = downsample(im, scale);
 if nargin<2
    im =double( imread('lena256.pgm'));
    scale = 2;
 end
 
[a b]   =       size(im);
new_r   =       1:scale:a;
new_c   =       1:scale:b;
im_new  =        im(new_r, new_c);
end