%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [baseLayer] = copy_and_downsample (src, d)
    d = 2^d;
    [w, h] = size(src);
    new_w  = 1:d:w;
    new_h  = 1:d:h;
    baseLayer = src(new_w, new_h);
end