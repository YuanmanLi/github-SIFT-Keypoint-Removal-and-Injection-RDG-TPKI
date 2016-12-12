%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I_new,injected_key, num_giveup, num_refine_out] = inject_keypoint(I_groundtruth,I,locations,injected_key, par)
%dbstop if all error
num_success      =   0;
num_giveup       =   0;
num_refine_out   =   0;
num_nofeasible   =   0; 
dis              =   3;
dog_rank         =   28 - 2*par.rank; %-1*(par.rank -1) + (27-par.rank)


max_stop_error  = par.stop_max_error;
cur_s                 = par.cur_s;
if nargin<1
    fprintf('no image input!\n')
end
if  size(locations) < 1
    fprintf('no injection  location\n');
    pause;
    return;
end
win_var    =     par.varwin;                                 %number of variable centered in win1.                                          %only consider the win1 block centered at curr_location. must bigger than w+4                                      
blkwin     =      [1 1];                              %inject keypoint at a specific location
win1       =   [floor(win_var(1)/2)+1, floor(win_var(2)/2)+1 ];                 
[num_points, ~] = size(locations);
I_new           = I;
[a b]           = size(I);
sum_err_sqr     = 0;
step            = 2^(par.cur_o-1);

for z =1:num_points
    giveup = 0;
    %%%%check if location is on the boundary%%%%%%%%%%%%%%%%%%%%
    current_point   = locations(z,1:2);
    space_boundary1 = floor(win1(1)/step)+1;
    space_boundary2 = floor(win1(2)/step)+1;
    if (current_point(1)+space_boundary1 > floor((a-1)/step)+1 | current_point(1)-space_boundary1 < 1 | current_point(2)+space_boundary2>floor((b-1)/step)+1 | current_point(2)-space_boundary2 < 1)
        fprintf('avoid the boundary keypoints\n');
        continue;
    end
    
    %%%%%%%%%%%%%%%%clip the large images for faster speed%%%%%%%%%%%%%%%%%%%%
    if par.cur_o ==1 || par.cur_o ==0
        clip_size = 32;
    else
        [clip_size,~] = size(I_groundtruth);
    end  
    cur_inj_loc = current_point;  
    temp1 = current_point(1);   temp2 = current_point(2);
    temp1 = ceil(temp1*(2^(par.cur_o-1))); temp2 = ceil(temp2*(2^(par.cur_o-1)));  %map the keypoint location to orginal image 
 
    %%%%check if there are any injected keypoint in the neighborhood%%%%%%%%
    if (IsNeighbor([temp1, temp2], injected_key,dis)>0)
       % fprintf('there are already some injected keypoints, give up\n');
        continue;
    end
    %%%%%%%%%%%end check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clip_indx_r =  max(1, temp1 - clip_size): min(a, temp1+clip_size);  %the position of clip image in original image
    clip_indx_c =  max(1, temp2 - clip_size): min(b, temp2+clip_size);
    clip_image  =  I_new(clip_indx_r, clip_indx_c);
    I = clip_image;
    step = 2^(par.cur_o-1);
    if par.cur_o ==0                                                                    %the upleft location of clip image in the upsmapled image
        clip_to_resize1    = clip_indx_r(1)/step-1;                          %not the interpolated pixel shold be in the even location
        clip_to_resize2    = clip_indx_c(1)/step-1;                          %so the groundtruth should be in the odd location
    else                                                                                        %the upleft location of clip image in the downsampled image
        clip_to_resize1    = ceil(clip_indx_r(1)/step);                   
        clip_to_resize2    = ceil(clip_indx_c(1)/step);                     
    end
    current_point(1)      = current_point(1) - clip_to_resize1+1;                  %location in the upsampled(downsampled) image
    current_point(2)      = current_point(2) - clip_to_resize2+1;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%calculate location of variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if par.cur_o ==0 || par.cur_o ==1    
         var_indx_r_ori = temp1 - floor(win_var(1)/2) : temp1 + floor(win_var(1)/2);                        %the variable location is the original image
         var_indx_c_ori = temp2 - floor(win_var(2)/2) : temp2 + floor(win_var(2)/2);      
    else
        var_indx_r1 = current_point(1) - floor(win_var(1)/2) : current_point(1) + floor(win_var(1)/2);
        var_indx_c1 = current_point(2) - floor(win_var(2)/2) : current_point(2) + floor(win_var(2)/2);       
        var_indx_r_ori = var_indx_r1 * step -1;                          %remap the variable location to original image
        var_indx_c_ori = var_indx_c1 * step-1;
    end
        var_indx_r  = var_indx_r_ori  - clip_indx_r(1) + 1;          %remap the variables location to clipped image
        var_indx_c  = var_indx_c_ori - clip_indx_c(1) + 1;
    %%%%%%%%%%%%%%construct R for 27 pixels%%%%%%%%%%%%%%%
    k = 0;                                                       % how many element in patch
    R = cell(1,blkwin(1)*blkwin(2));                              %9 matrix extractor for a given 3*3 patch
    cur_o_size =  size(I)./(2^(par.cur_o-1));
    for i = 1:blkwin(2)
        for j = 1: blkwin(1)
            k = k +1; 
            c_loc = [(current_point(1)-floor(blkwin(2)/2) + i -1),  (current_point(2)-floor(blkwin(1)/2) + j -1)]; %get all the locations in the patch
            R{k} = getR(c_loc, [3 3],cur_o_size);                                                                                                                    %get the patch 3*3 patch extracting matrix
        end
    end
    %%%%%%%%%%parameter setting %%%%%%%%%%%%
    f.nlevels = 3;
    f.octave = 1;
    f.o_min = par.o_min;
    f.s_min = -1;

    f.sigman = 0.5;
    f.sigmak = 2.0^(1.0/f.nlevels);
    f.sigma0 = 1.6 * f.sigmak;
    f.dsigma0 = f.sigma0 * sqrt(1.0 - 1.0/(f.sigmak*f.sigmak));

    f.peak_thresh = 0;                                                                      %this is used for finding local minma or maxma
    f.edge_thresh= 10;
    f.norm_thresh=0;
    f.magnif = 3;
    f.windowSize=2;
    
    if f.o_min==-1
        %fprintf('handle upsample\n');
        I_resize =  copy_and_upsample (I, 1);
    elseif f.o_min>0
        fprintf( 'handle downsample');
    else
        I_resize = I;                                                                                    %just the original image
    end
    sa  = f.sigma0* (f.sigmak^ f.s_min);
    sb = f.sigman * 2.0^(-f.o_min);

%%handle the first level 
if(sa > sb)
    sd = sqrt(sa*sa - sb*sb);
end
G_cell = cell(1,6);
Sig_cell = zeros(1,6);
Sig_cell(1) = sd;
[x]=getGuassianSize(Sig_cell(1));    
G_cur = fspecial('gaussian', [x x], Sig_cell(1)); 
G_cell{1} = G_cur;
idx = 1;
for  s = f.s_min+1 : f.nlevels - f.s_min
    idx = idx+1;
    sd = f.dsigma0 * (f.sigmak^ s) ;
    Sig_cell(idx) = sd;
    [x]=getGuassianSize(Sig_cell(idx));    
    G_cur = fspecial('gaussian', [x x], Sig_cell(idx)); 
    G_cell{idx}  = G_cur;
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
       G2_I =  imfilter(G1_I,G_cell{2},'replicate');
       G3_I =  imfilter(G2_I,G_cell{3},'replicate');
       G4_I =  imfilter(G3_I,G_cell{4},'replicate');
     end 
    G5_I   =  imfilter(G4_I,G_cell{5},'replicate');
    G6_I   =  imfilter(G5_I,G_cell{6},'replicate');
   
    D1 =  G2_I-G1_I;
    D2 =  G3_I-G2_I;
    D3 =  G4_I-G3_I;
    D4 =  G5_I-G4_I; 
    D5 =  G6_I-G5_I;    
    D_cell = {D1(:), D2(:), D3(:), D4(:), D5(:)};
    Ks     = zeros(27, blkwin(1)*blkwin(2), s); 
    cent   = floor(blkwin(1)*blkwin(2)/2)+1;
    key_Max= 0;                                                       %if key_Max =0, give up
    
for L = cur_s: cur_s 
        K = zeros(27, blkwin(1)*blkwin(2));                %every column is the neighbor of current location 
        k = 0;                                                                  % how many element in patch
        for i = 1:blkwin(1)
            for j = 1: blkwin(2)
                k = k +1;
                K(:,k) = [R{k}*D_cell{L}; R{k}*D_cell{L+1}; R{k}*D_cell{L+2}];      %D1, D2, D3 are the DOGs
                Ks(:, k, L) =   [R{k}*D_cell{L}; R{k}*D_cell{L+1}; R{k}*D_cell{L+2}];
            end
        end
        %%%%%%check keypoint%%%%%%%%%%%%%%%%%%%
        if L ==cur_s
            patch = Ks(:, cent, L);                         %for injection, cent and L should one
            cur_Dog = patch(14);
            dif_patch = sign( cur_Dog - patch );  
            if sum(dif_patch) <=-dog_rank                      %%%a good location should be at the end(dog vaule should be eigher small or big)
                key_Max = -1;
            elseif sum(dif_patch) >= dog_rank
                key_Max = 1;
            else 
                giveup = 1;
            end
        end 
end
    if giveup ==1
         continue;
    end
  
  %%%%%%%%%begin optimization%%%%%%%%%%
  I_temp = I_groundtruth(var_indx_r_ori, var_indx_c_ori);
  I_col  = double( I_temp(:));
  
  opts = optimset('LargeScale','off','display','off','Algorithm','active-set','TolCon',1e-8,'MaxFunEvals',10000);
  [x,fval,exitflag]=fmincon('obj_fun_inject',I_col,[],[],[],[],[],[],'con_fun_inject',opts,I_col, I,G_cell,R,blkwin, win_var,var_indx_r, var_indx_c, key_Max,par);

  if exitflag == -2
     num_nofeasible = num_nofeasible + 1;
  end
  %%%%%%%%%% update the win1 patch%%%%%%%%%%%%%%%%%
  curr_err_sqr = sum((x-I_col)'*(x-I_col));
  %%%%%%%give up when error is too big%%%%%%%%%%%%%%%
  if max(abs(x-I_col)) < max_stop_error              
        sum_err_sqr = sum_err_sqr+curr_err_sqr;
        fprintf('curr_err_sqr = %.2f, max_err = %.2f, sum_err_sqr=%.2f\n', curr_err_sqr, max(abs(x-I_col)),sum_err_sqr);
        dif = x - I_col;
        x(dif < 0) = floor(x(dif < 0));
        x(dif > 0) = ceil(x(dif > 0));
        x=reshape(x, win_var);                                   %updated 
        x_int = uint8(x);
    %%%%%%%%%%%%check whether this keypoint has been successfully injected%%%%
    I_temp = I_new;
    I_temp(var_indx_r_ori, var_indx_c_ori) = x_int;   %for clip image
    flag = checkInjection(cur_inj_loc, I_temp,par);
    if(flag == 1)
        fprintf('injecting success %d, %d\n',temp1, temp2);
        num_success = num_success + 1;
        I_new = I_temp;
        injected_key = [injected_key; [temp1, temp2]];
        fprintf('%d of %d \n ', z, num_points);
    else
        num_refine_out = num_refine_out +  1;
        fprintf('inject fail, refinement filtering out,%d of %d \n',z, num_points);
        continue;
    end
    %%%%%%%%%%%end check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  else
      num_giveup = num_giveup + 1;
      fprintf('error is big, give up\n')
      fprintf('curr_err_sqr = %.2f, max_err = %.2f\n', curr_err_sqr, max(abs(x-I_col)));
 end
end
%%%%%%%%%%%check the remained injected keypoints, and exclude keypoints
%%%%%%%%%%%dismissed%%%%%%%%%%%%%%%%%%%%%
injected_key = refine_injec_locs(injected_key, I_new,par);
fprintf('stop here\n');
fprintf('-----------------inject %d keypoints \n', num_success );
fprintf('-----------------totally injected %d keypoints \n', size(injected_key,1) );
fprintf('-----------------totally giveup %d error big injection \n', num_giveup );
fprintf('  ----------------totally given up  %d refine_out injection\n', num_refine_out);
fprintf('  ----------------toally  found  %d not feasible keypoints\n', num_nofeasible);
end

function flag = IsNeighbor(point, neighbor,dis)
    flag = 0;
    difx = abs(neighbor(:,1)-point(1));
    dify = abs(neighbor(:,2)-point(2));
    indx = difx<=dis & dify<=dis;
    if(sum(indx)>0)
        flag = 1;
    end
end

function flag = checkInjection(injec_loc, I_inject,par)
    flag = 0;
    dis  = 1;                            %because C and matlab have a small margin of error £¨set dis=0 just ignore such phenomenon£©
    [~,~,keypoints_src,~]       =           vl_sift(single(I_inject),'Octaves',par.octaves,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',par.cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
    keypoints_src =  double(keypoints_src);
    difx = abs(keypoints_src(1,:) - injec_loc(1));
    dify = abs(keypoints_src(2,:) -injec_loc(2));
    ind = difx<=dis & dify<=dis;
    if sum(ind)>0
        flag = 1;
    end
end   