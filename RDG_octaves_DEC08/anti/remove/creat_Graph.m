%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m_adj] =  creat_Graph(data, S, win_w, key_Max)

m_indx   =   1: win_w*win_w;
m_indx   =    reshape(m_indx, win_w, win_w);
len_win  =    win_w *win_w;
t_indx(:,:,1)  =   m_indx';
for i = 1:S+1
    t_indx(:,:,i+1)  = t_indx(:,:,i)  + len_win;                % the index of every Dog Layer
end
i       =        2 : 4; 
j       =        2 : 4;
s       =        2 : S+1;
b1      =        win_w-2;
b2      =        win_w-2;

G_indx  =  zeros((win_w-2) ^2 * S,27);                                 % win_w-2 means the center patch in which there is no keypoint. 
k = 1;
for L  = 1:3 
    for k1 = 1 : b1
        for k2 = 1 : b2
                indx       =   t_indx (i+k1-2, j+k2-2, s+L-2);    %  be careful. even the block 27*27 is scanned as row, but the 27 adjacent pixels of each is scanned as columns. 
                indx       =   indx(:);
                G_indx(k,:)=   indx';
                k          =   k + 1;
        end
    end
end

max_dis     =    10^6; 
m_adj       =     max_dis * ones(win_w * win_w * (S+2)); 
pos_key     =     G_indx(:,14);                                                        % the position of center pixels
%%for test %%
[len,~]    =  size(G_indx);
for i = 1:len
       s            =    floor((i-1)/(b1*b2))+1;                                                    % because the indx in matlab begins from 1. so i -1
       j            =    i-(s-1)*(b1*b2);
       temp     =    data(:,j,s);
       temp     = temp - temp(14);                                                    % the order means the distance from the center to neighbor. 
       temp(14) = max_dis;                                                              %because no distance to itself
       m_adj(pos_key(i), G_indx(i,:)) = temp;
end


