%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c , ceq]  =  con_fun_next_o(x, I_col,I,  G_cell, R,win,  s, win_var,var_indx_r, var_indx_c, max_min_loc, par)

    y = reshape(x, win_var);
    I(var_indx_r, var_indx_c) = y;
    
    if par.o_min == -1
        I_resize =  copy_and_upsample (I, 1);
        I        =  I_resize;
    end
    G1_I = imfilter(I ,G_cell{1},'replicate');
    G2_I =  imfilter(G1_I,G_cell{2},'replicate');
    G3_I =  imfilter(G2_I,G_cell{3},'replicate');
    G4_I =  imfilter(G3_I,G_cell{4},'replicate');

     for i = par.o_min+2:par.cur_o                                                     %for higher octave
       G1_I = copy_and_downsample (G4_I, 1);         %changed to the best layer
       G2_I =  imfilter(G1_I,G_cell{2},'replicate');
       G3_I =  imfilter(G2_I,G_cell{3},'replicate');
       G4_I =  imfilter(G3_I,G_cell{4},'replicate');
     end 
    G5_I =  imfilter(G4_I,G_cell{5},'replicate');
    G6_I =  imfilter(G5_I,G_cell{6},'replicate');
   
    D1 =  G2_I-G1_I;
    D2 =  G3_I-G2_I;
    D3 =  G4_I-G3_I;
    D4 =  G5_I-G4_I;
    D5 =  G6_I-G5_I;
    
    D_cell = {D1(:), D2(:), D3(:), D4(:), D5(:)};
    b1     =    win(1); b2 = win(2);
    blk_size = b1*b2;
    %%%%%%%%%%%%%%construct K%%%%%%%%%%%%%%%
    %% every column in K denotes the 9 neighborhood of the corresponding pixels in win. 
   if ismember(par.cur_o, par.o_path)     
        thr = par.con_threthold;                                                         %to aviod equal constraint      
   else
       thr =0;
   end

   c       = zeros(blk_size*s*2,1);
   for L    = 1: s 
        K     =    zeros(27, win(1)*win(2));                %every column is the neighbor of current location 
        k     =    0;                                                                  % how many element in patch
        for i = 1:b2
            for j = 1: b1
                k      = k +1;
                K(:,k) = [R{k}*D_cell{L}; R{k}*D_cell{L+1}; R{k}*D_cell{L+2}];      %D1, D2, D3 are the DOGs
            end
        end

        for i = 1:b1*b2
            temp_loc        =       max_min_loc((L-1)*blk_size+i, :);
            max_loc         =       temp_loc(1); 
            min_loc         =        temp_loc(2);
            curr_K          =        K(:,i);

            c((L-1)*blk_size*2+2*i-1)     =  curr_K(14) - curr_K(max_loc) + thr;        % max>= x >= min
            c((L-1)*blk_size*2+2*i)       =  curr_K(min_loc) - curr_K(14) + thr;
        end
    
   end
    ceq =[];
end
