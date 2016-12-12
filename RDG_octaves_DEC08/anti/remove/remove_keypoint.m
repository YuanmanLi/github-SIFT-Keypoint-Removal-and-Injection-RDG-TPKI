%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I_new, bd_keypoint, num_giveup, num_goodpath, num_abnormal] = remove_keypoint(I_groundtruth,I,locations, par)
num_giveup     = 0;
num_goodpath   = 0;
num_abnormal   = 0;
num_nofeasible = 0; 
max_stop_error = par.stop_max_error;
if nargin<1
    error('no image input!\n')
end
if  size(locations) < 1
    fprintf('no keypoints location\n');
    pause;
    return;
end                                        
[num_points, ~] =  size(locations);
I_new       = I;
[a, b]      = size(I);
sum_err_sqr = 0;
step        = 2^(par.cur_o-1);  %resize the input image for different octaves


bd_keypoint  = [];
for z =1:num_points
    
    blkwin  =  par.blkwin;                              
    win1    =  [floor(blkwin(1)/2)+1, floor(blkwin(2)/2)+1 ];                                           %only consider the win1 block centered at curr_location. must bigger than w+4   
    win_var =  par.varwin;      %number of variable centered in win1.  
 
    current_point = locations(z,1:2);
    cur_s         =  locations(z,3);  
    fprintf('%d of %d (location: [%d,%d,%d])\n ', z, num_points,current_point(1),current_point(2),cur_s);
    space_boundary1 = win1(1);
    space_boundary2 = win1(1);
    %%%%%%%%%%keypoint location%%%%%%%%%%%%%
    if (current_point(1)+space_boundary1 > floor((a-1)/step)+1 | current_point(1)-space_boundary1 < 1 | current_point(2)+space_boundary2>floor((b-1)/step)+1 | current_point(2)-space_boundary2 < 1)
         fprintf('handle the boundary keypoints\n');
         bd_keypoint = [bd_keypoint; current_point];
         blkwin      =     [1 1];                             
         win_var     =   [3 3];
         if par.cur_o==0 && (current_point(1) <3 ||current_point(2) <3 ||current_point(1)>floor((a-1)/step)-1 || current_point(2)>floor((b-1)/step)-1) 
            fprintf('cannot handle the boundary keypoints\n');
            continue;
         end
         if par.cur_o ==1 && (current_point(1)==1 || current_point(2)==1)
             fprintf('abnomal! cannot handle the boundary keypoints\n');
             continue;
         end
    end
    
    %%%%%%%%%%%%%%%%clip the large images for faster speed%%%%%%%%%%%%%%%%%%%%
    if par.cur_o ==1 || par.cur_o ==0
        clip_size = 16; %change it as 32 for accuracy
    else
        clip_size = max(size(I_groundtruth));
    end  
    temp1 = current_point(1); temp2 = current_point(2);
    temp1 = ceil(temp1*(2^(par.cur_o-1))); temp2 = ceil(temp2*(2^(par.cur_o-1)));  %map the keypoint location to orginal image 
 
    clip_indx_r  =  max(1, temp1 - clip_size): min(a, temp1+clip_size);            %the position of clip image in original image
    clip_indx_c  =  max(1, temp2 - clip_size): min(b, temp2+clip_size);
    clip_image   =  I_new(clip_indx_r, clip_indx_c);
    I            =  clip_image;

    step = 2^(par.cur_o-1);
    if par.cur_o ==0                                                         
        clip_to_resize1    = clip_indx_r(1)/step-1;                          
        clip_to_resize2    = clip_indx_c(1)/step-1;                         
    else                                                                     
        clip_to_resize1    = ceil(clip_indx_r(1)/step);                   
        clip_to_resize2    = ceil(clip_indx_c(1)/step);                     
    end
    current_point(1)      = current_point(1) - clip_to_resize1+1;                  %location in the upsampled(downsampled) image
    current_point(2)      = current_point(2) - clip_to_resize2+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % I = I_new;
    %%%%%%%%%%calculate location of variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if par.cur_o ==0 || par.cur_o ==1    
         var_indx_r_ori = temp1 - floor(win_var(1)/2) : temp1 + floor(win_var(1)/2);                        %the variable location is the original image
         var_indx_c_ori = temp2 - floor(win_var(2)/2) : temp2 + floor(win_var(2)/2);      
    else
        var_indx_r1 = current_point(1) - floor(win_var(1)/2) : current_point(1) + floor(win_var(1)/2);
        var_indx_c1 = current_point(2) - floor(win_var(2)/2) : current_point(2) + floor(win_var(2)/2);       
        var_indx_r_ori = var_indx_r1 * step -1;   %remap the variable location to original image
        var_indx_c_ori = var_indx_c1 * step-1;
    end
        var_indx_r  = var_indx_r_ori  - clip_indx_r(1) + 1;              
        var_indx_c  = var_indx_c_ori - clip_indx_c(1) + 1;
    
    %%%%%%%%%%%%%%construct R for 27 pixels%%%%%%%%%%%%%%%
    k = 0;                                                                                              % how many element in patch
    R = cell(1,blkwin(1)*blkwin(2));                                                                                   %9 matrix extractor for a given 3*3 patch
    cur_o_size =  size(I)./(2^(par.cur_o-1));
    for i = 1:blkwin(2)
        for j = 1: blkwin(1)
            k = k +1;
            c_loc = [(current_point(1)-floor(blkwin(2)/2) + i -1),  (current_point(2)-floor(blkwin(1)/2) + j -1)]; %get all the locations in the patch
            R{k} = getR(c_loc, [3 3],cur_o_size);                                                                                                                    %get the patch 3*3 patch extracting matrix
        end
    end
    %%%%%%%%%%parameter setting for VLSIFT %%%%%%%%%%%%
    f.nlevels = 3;
    f.octave  = 1;
    f.o_min   = par.o_min;
    f.s_min   = -1;
    f.sigman  = 0.5;
    f.sigmak  = 2.0^(1.0/f.nlevels);
    f.sigma0  = 1.6 * f.sigmak;
    f.dsigma0 = f.sigma0 * sqrt(1.0 - 1.0/(f.sigmak*f.sigmak));
    f.peak_thresh = 0;                                                                      %this is used for finding local minma or maxma
    f.edge_thresh = 10;
    f.norm_thresh = 0;
    f.magnif      = 3;
    f.windowSize  = 2; 
    if f.o_min==-1
        I_resize =  copy_and_upsample(I, 1);
    elseif f.o_min>0
        fprintf( 'handle downsample');
    else
        I_resize = I;                                                                                    %just the original image
    end
    sa  = f.sigma0 * (f.sigmak^ f.s_min);
    sb  =  f.sigman * 2.0^(-f.o_min);
    %%handle the first level 
    if(sa > sb)
        sd = sqrt(sa*sa - sb*sb);
    end
    guassian.sigma = sd;

    G_cell   = cell(1,6);
    Sig_cell = zeros(1,6);
    Sig_cell(1) = sd;
    [x]      = getGuassianSize(Sig_cell(1));    
    G_cur    = fspecial('gaussian', [x x], Sig_cell(1)); 
    G_cell{1}= G_cur;
    idx = 1;
    for  s  = f.s_min+1 : f.nlevels - f.s_min
        idx = idx+1;
        sd  = f.dsigma0 * (f.sigmak^ s) ;
        Sig_cell(idx) = sd;
        [x] = getGuassianSize(Sig_cell(idx));    
        G_cur = fspecial('gaussian', [x x], Sig_cell(idx)); 
        G_cell{idx}= G_cur;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%return the location of maximum and minimum value%%%%%%%%%%%%%%
    s    =  f.nlevels;
    G1_I = imfilter(I_resize,G_cell{1},'replicate');
    G2_I =  imfilter(G1_I,G_cell{2},'replicate');
    G3_I =  imfilter(G2_I,G_cell{3},'replicate');
    G4_I =  imfilter(G3_I,G_cell{4},'replicate');

     for i = par.o_min+2:par.cur_o                                                     %for higher octave
       G1_I = copy_and_downsample (G4_I, 1);         %changed to the best layer
%        sa = f.sigma0 * (f.sigmak^ f.s_min) ;
%        sb = f.sigma0 *(f.sigmak^( s_best - f.nlevels)) ;
%        if sa > sb
%          G1_I =  imfilter(G1_I,G_cell{1},'replicate');
%        end 
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
    max_min_loc = zeros(s*blkwin(1)*blkwin(2),2);    %max_min_loc = (3*5*5)*2
    max_min_loc_s = zeros(s, s*blkwin(1)*blkwin(2),2);
    Ks  = zeros(27, blkwin(1)*blkwin(2), s); 
    cent = floor(blkwin(1)*blkwin(2)/2)+1;
    key_Max = 0;                                      %if key_Max =0, then it's not the keypoint, if =-1, it's the minimumn, if =1, it's the maximumn
for L = 1: s 
        K = zeros(27, blkwin(1)*blkwin(2));           %every column is the neighbor of current location 
        k = 0;                                        % how many element in patch
        for i = 1:blkwin(1)
            for j = 1: blkwin(2)
                k = k +1;
                K(:,k)      =   [R{k}*D_cell{L}; R{k}*D_cell{L+1}; R{k}*D_cell{L+2}];      %D1, D2, D3 are the DOGs
                Ks(:, k, L) =   [R{k}*D_cell{L}; R{k}*D_cell{L+1}; R{k}*D_cell{L+2}];
            end
        end
        %%%%%%check keypoint%%%%%%%%%%%%%%%%%%%
        if L ==cur_s
            patch = Ks(:, cent, L);
            cur_Dog = patch(14);
            dif_patch = sign( cur_Dog - patch );  
            if sum(dif_patch) == -26
                key_Max = -1;
            elseif sum(dif_patch) == 26
                key_Max = 1;
            else 
                disp('it is not the keypoint');
            end
        end
end

%%%%%%%%%create Graph %%%%%%%%%%%%%
  if key_Max ~= 0     %current is the keypoint
          block_size   =  (blkwin(1)+2) * (blkwin(2)+2);
          start_indx   =   block_size * cur_s + ceil(block_size/2);
          m_adj        =   creat_Graph(Ks,s, blkwin(1)+2, key_Max);
          [path, ~ ,~, ~, b_goodpath]    =   shortpath(m_adj, start_indx, key_Max);
          num_goodpath =    num_goodpath + b_goodpath;
          [p_map]      =    path_map(path, blkwin(1)+2);
          [path_len, ~]= size(p_map);
          if path_len <  2 && ~((cur_s == 1&& p_map(3)<=9) ||  (cur_s == s && p_map(3)>=19))                 %no keypont will be generated at the lowest layer and highest layer
               disp('abnormal~!\n');
          end
  else 
          disp('it is not keypoint!\n');
          continue;   %still try to remove, because there is a little difference to vl_feat_c
  end
for L = 1: s 
    %%%%%%%%%%%%%do not use path algorithm%%%%%%%%%%%%%%%%%%%%%%%  
    if ismember(par.cur_o, par.o_path)            
        s_path   =   p_map(p_map(:,1)==L,:);
        K            =    Ks(:, :, L);
        %%write as two functions
        min_loc = 0; max_loc = 0;
        for i = 1 : blkwin(1)*blkwin(2) 
            curr_Dog       =   K(:,i);                                   
            [~, sort_indx] =   sort(curr_Dog); 
            min_loc        =   sort_indx(1);
            max_loc        =   sort_indx(end);
            if ~isempty(s_path) 
               next_p     =   s_path(:,3);
               temp_indx  =   find(i == s_path(:,2)-blkwin(1)*blkwin(2) );
               temp_extre =   next_p(temp_indx);
               if  key_Max  == -1 && ~isempty(temp_extre) 
                   min_loc = temp_extre;
                   s_path(temp_indx,:) = [];                                         
               elseif key_Max == 1 && ~isempty(temp_extre) 
                   max_loc = temp_extre;
                   s_path(temp_indx,:) = [];                                         
               end
            end       
        if max_loc == 14 || min_loc ==14
            num_abnormal = num_abnormal + 1;
            disp('position of maximum and minimumn can not be 14\n')       %perhaps 1 or 3 layers both have keypoint 
        end      
        max_min_loc((L-1)*blkwin(1)*blkwin(2) + i, :) = [max_loc, min_loc]; 
        end
        if ~isempty(s_path)
            error('there still some path(s) not be used\n');
        end
      %%%%%%%%%%%%%do not use path algorithm%%%%%%%%%%%%%%%%%%%%%%%  
    else                                     
        K            =    Ks(:, :, L);
        %%write as two functions
        min_loc = 0; max_loc = 0;
        for i = 1 : blkwin(1)*blkwin(2)
            curr_Dog = K(:,i);                                   
            [~, sort_indx] = sort(curr_Dog); 
            min_loc          =       sort_indx(1);
            max_loc         =       sort_indx(end);
            if sort_indx(1) == 14                           %equal to 14 mean current location is keypoint
                min_loc       =   sort_indx(2);
            end
            if sort_indx(end) == 14
                max_loc      =   sort_indx(end-1);
            end
            max_min_loc((L-1)*blkwin(1)*blkwin(2) + i, :) = [max_loc, min_loc]; 
        end
    end
end
  
  
  %%%%%%%%%begin optimization%%%%%%%%%%
  I_temp = I_groundtruth(var_indx_r_ori, var_indx_c_ori);
  I_col=double( I_temp(:));
  opts = optimset('LargeScale','off','display','off','Algorithm','active-set','TolCon',1e-8,'MaxFunEvals',500,'MaxIter',50);
  [x,fval,exitflag]=fmincon('obj_fun',I_col,[],[],[],[],[],[],'con_fun_next_o',opts,I_col, I,G_cell,R, blkwin,  f.nlevels, win_var,var_indx_r, var_indx_c, max_min_loc, par);
  
  if exitflag == -2
     num_nofeasible = num_nofeasible + 1;
  end

  %%%%%%%%%% update the win1 patch%%%%%%%%%%%%%%%%%
  curr_err_sqr = sum((x-I_col)'*(x-I_col));
  %%%%%%%give up when error is too big%%%%%%%%%%%%%%%
  if max(abs(x-I_col)) < max_stop_error              
        sum_err_sqr = sum_err_sqr+curr_err_sqr;
        fprintf('curr_err_sqr = %.2f, max_err = %.2f, sum_err_sqr=%.2f\n', curr_err_sqr, max(abs(x-I_col)),sum_err_sqr);
        x    =  reshape(x, win_var);                       %updated 
        I_new(var_indx_r_ori, var_indx_c_ori) = x;   %for clip image
  else
      num_giveup = num_giveup + 1;
      fprintf('error is big, give up\n')
      fprintf('curr_err_sqr = %.2f, max_err = %.2f\n', curr_err_sqr, max(abs(x-I_col)));
 end
end
fprintf('stop here\n');
fprintf('-----------------totally giveup %d removed keypoint \n', num_giveup );
fprintf('  ----------------totally found  %d good path\n', num_goodpath);
fprintf('  ----------------totally found  %d abnormal pixels ( perhaps 1 or 3 layers both have keypoint )\n', num_abnormal);
fprintf('  ----------------toally  found  %d not feasible keypoints\n', num_nofeasible);
% a;
