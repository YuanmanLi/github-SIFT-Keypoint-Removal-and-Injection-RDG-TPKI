%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R] = getR(curr_location, win, wh)
w = wh(1); h = wh(2);
n          =  w*h;
pic_idx    = 1:n;
pic_inx_m  = reshape(pic_idx, w, h);
b1 = win(1); b2 = win(2);

r = curr_location(1)-floor(b1/2)    :   curr_location(1)+floor(b1/2);
c =  curr_location(2)-floor(b2/2)    :   curr_location(2)+floor(b2/2);
patch_idx_m = pic_inx_m(r, c);

R = sparse(1:b1*b2, patch_idx_m(:), ones(1, b1*b2), b1*b2, n, b1*b2 );

