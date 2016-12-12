%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [baseLayer] = copy_and_upsample (src, d)
    d            =   2^d;
    [w, h]       =    size(src);
    interp       =  (src(1:w, 1:h-1)+src(1:w,2:h))/2;  % interp: w*(h-1)
    new_w   = 1:d:2*w;
    new_h   = 1:d:2*h;
    baseLayer                       =   zeros(2*w, 2*h);
    baseLayer(new_w, new_h)         =   src;
    baseLayer(new_w, new_h(1:end-1)+1) =interp;
    baseLayer(:, 2*h) =    baseLayer(:, 2*h-1);                           %the boundary
    
    interp1     =  (baseLayer(new_w(1:end-1), :)+baseLayer(new_w(2:end),:))/2;  % interp: w*(h-1)
    baseLayer(new_w(1:end-1)+1, :) =    interp1;
    baseLayer(2*w,:) =    baseLayer(2*w-1,:) ;                           %the boundary
end