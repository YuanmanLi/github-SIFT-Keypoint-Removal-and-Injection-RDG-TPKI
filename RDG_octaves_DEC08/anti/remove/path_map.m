%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_map] = path_map(path, block_width)

	block_size          =       block_width^2;
    block_cen_size      =       (block_width -2) * (block_width -2);
    len                 =        length(path);
    p_map               =       zeros(len-1, 3);
    temp                =       1:block_size;
    temp                =       reshape(temp, block_width, block_width);
    indx_cen            =       temp(2:block_width-1, 2:block_width-1);               %every layer, 9pixel
    
    for i   =   1:len-1
              cen         =       path(i);
              cur_s       =        floor((cen-1)/block_size);
              M           =       (1: block_size) + cur_s*block_size;
              M           =       reshape(M, block_width, block_width);
              M = M';
              r           =       floor((cen -  cur_s*block_size -1) / block_width);
              c           =       cen -  cur_s*block_size - r * block_width;
              C           =       M(r :  r+2, c-1 : c+1);
              C_s         =        [C(:)-block_size; C(:); C(:)+block_size];
              indx_s1     =        find(C_s == path(i+1));
              indx_s2     =        find(indx_cen == mod(cen, block_size));
              p_map(i,:)  =        [cur_s,indx_s2 + block_cen_size, indx_s1];
    end
    
  
  

  
  
  
