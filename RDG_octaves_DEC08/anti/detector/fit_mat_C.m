%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%fit_mat_C is use to match the location of the result of c and matlab. 
%vl_feat used the float(single), while in matlab the single typle can't be
%used when doing the smooth and optimization
%scr: the location detected by vl_feat, 
%des: the location detected by matlab
function [keypoints] = fit_mat_C(src, des)
    src    = double(src);
    [m n]  = size(src);
    [m1 n1]= size(des);
    max_dis= 3;
    keypoints = 3*ones(m,n);
    for i = 1:m
        cur_key = src(i,1:2);
        dif1    = abs(des(:,1:2) - repmat(cur_key,m1,1));
        dif     = dif1(:,1)+dif1(:,2);
        [~,indx]= min(dif);
        keypoints(i,:) = des(indx(1),:);          
    end 
end