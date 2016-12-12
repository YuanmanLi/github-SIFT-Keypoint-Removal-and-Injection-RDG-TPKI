%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%% This is the implementation of the paper: "SIFT Keypoint Removal via Directed Graph Construction for Color Images"
%% The code also provides the component of injection algorithm TPKI. 
%% Use the code, please cite the following paper£º
%% Y. Li, J. Zhou, A. Cheng, X. Liu, and Y. Y. Tang, 
%% ¡°SIFT keypoint removal and injection via convex relaxation,¡± IEEE TIFS, vol. 11, no. 8, Aug 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function [I_cur,key_removed] =  main_remove(I,I_cur, injec_loc, utilPara,savepath, im_idx)
 %%%%%%paras%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.max_loo                  =    utilPara.max_loo ;                                                  %max loop                                                                              
par.con_threthold            =    0.01*10^(-0);     
par.o_path                   =    [0 1 2];             

par.save_im                  =    utilPara.save_im;                 %whether save image in each loop
par.zero_threshold           =    utilPara.zero_threshold ;    
par.peak_threshold           =    utilPara.peak_threshold;                                      %parameter for vl_feat constrast¡£ used in detect_keypoint.m
par.edge_threshold           =    utilPara.edge_threshold;
par.octaves                  =    utilPara.octaves;
par.o_min                    =    utilPara.o_min;
par.re_octaves               =    utilPara.re_octaves;

if par.o_min ==-1
 %   par.stop_max_error       =   15;   %set larger than 50 to speed up
    par.active_errorUpdate_iter= 10;   
    con_thre                 =   [0.01*10^(0),0.01*10^(0),0,0];                                  %for optimization. only work for 
    blkwins                  =   [5 5;5 5; 3 3; 3 3];
    varwins                  =   [3 3;7 7; 5 5; 5 5];                                           %number of variable centered in win1. 
    par.stop_max_errors      =   [15;15;15;15];                                                  %allowed maximum error
    max_error_updaterates    =   [2; 1.5; 1; 0.5];                                            %add max_error by max_error_updaterate from par.active_errorUpdate_iter
else
    par.active_errorUpdate_iter= 5;  
    con_thre                 =   [0.01*10^(-1)];     %for optimization. only work for 
    blkwins                  =   [5 5];              %change the size randomly might make the results better.                   
    varwins                  =   [7 7];              %smaller value for speed, higher value for quality                             
    par.stop_max_errors      =   [15];                                                  %allowed maximum error
    max_error_updaterates    =   [3];                %add max_error by max_error_updaterate from par.active_errorUpdate_iter
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loo                    =    1;
key_removed            =    zeros(par.max_loo, 7);
key_removed(loo, 1:2)  =    [10 1];   % nothing useful, just set the initial condition

[f_old,d_old] = vl_sift(single(I),'Octaves',par.octaves ,'Levels',3, 'FirstOctave',par.o_min,'ThreEqual',par.zero_threshold,'PeakThresh', 4, 'EdgeThresh', par.edge_threshold);

while  loo <= par.max_loo   
       fprintf('\ncurrent loo=%d, (imdx=%d)\n',loo,im_idx);
       %%%%%%%%%%%%iteration setting %%%%%%%%%%%%%%%
       con_thre        =       con_thre * 0.95;  %as the left key points are hard to remove (max error will be large)
       if loo   >=   par.active_errorUpdate_iter && loo<=50
            par.stop_max_errors = par.stop_max_errors + max_error_updaterates; 
       end
       %%%%%%%%%%%%end setting %%%%%%%%%%%%%%%%%%
   for o_loop = par.re_octaves(1):par.re_octaves(2)
       cur_octave          =    o_loop-par.o_min;
       par.con_threthold   =    con_thre(cur_octave);     
       par.blkwin          =    blkwins(cur_octave,:);
       par.varwin          =    varwins(cur_octave,:);                                          %number of variable centered in win1. 
       par.stop_max_error  =    par.stop_max_errors(cur_octave);
       new_o_loop          =    o_loop;
       if loo <= par.max_loo*0.7 && o_loop > 1                                  
            continue;
       elseif loo<=par.max_loo*0.8 && o_loop >2                                   
            continue;
       elseif loo>par.max_loo*0.7
           par.o_path                  =       [0 1];                            
       elseif loo>par.max_loo*0.8
           par.o_path                  =       [0];                                
       end
       if loo > par.max_loo-4 &&o_loop~=par.re_octaves(1) 
           continue; 
       elseif loo>par.max_loo -7 && loo <= par.max_loo-4 && o_loop~=1  
           continue;
       elseif loo>par.max_loo -9 && loo <= par.max_loo-7 && o_loop ~= 2; 
           continue; 
       elseif loo > par.max_loo-4
           if (rand(1)>0.5) %change the size randomly
              %par.blkwin =blkwins(cur_octave,:)-2;
              par.varwin =varwins(cur_octave,:)-2;     
           end
       end
       par.cur_o                       =   new_o_loop;
       [f_new,d_new,keypoints,keyrep]  =   vl_sift(single(I_cur),'Octaves',par.octaves ,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',par.cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
       keypoints    =   keypoints';
%       [~,idx]      =   sort(keypoints(:,1));   %should not sort, because here use keyrep                       
%       keypoints    =   keypoints(idx,:);  
       [keypoints1] =   detect_keypoint_nextO(I_cur, par);
       keypoints    =   fit_mat_C(keypoints, keypoints1);
%        fprintf('before matching refine: %d\n', size(keypoints,1));
        
       [matches,~]  =   vl_ubcmatch(d_old, d_new) ;
       [r_Matches]  =   JudgeMat(f_old, f_new, matches);  
       [matched_keypoints]  =  calc_matched_keypoint(keyrep, r_Matches, par.octaves);
       new_key_indx         =  unique(matched_keypoints{cur_octave});
       keypoints            =  keypoints(new_key_indx,:);
%       fprintf('after matching refine: %d\n', size(keypoints,1));
        
        %%%%%%%%%%%%%%%exclude injected ones%%%%%%%%%%%
       [keypoints]  =  refine_remove_keys(keypoints, injec_loc);
%        fprintf('after remove-injection refine: %d\n', size(keypoints,1));
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        key_removed(loo, 1)  =  size(keypoints,1);
%        key_removed(loo, 6)  =  size(r_Matches,2);
        %%%%%%%%%%%%%%%%removing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       [I_new bd_keypoint]  =  remove_keypoint(I_cur,I_cur,keypoints, par); 
         %%%%%%%%%%%%%%%%end removing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if par. save_im == 1 
                save_path             = strcat(savepath,int2str(im_idx),'_iter', int2str(loo), '.png');
                fprintf('save image here\n');
                I_new_save            = I_new;
                I_dif                 = I_new - I;
                I_new_save(I_dif < 0) = floor(I_new_save(I_dif < 0));%use random might be better
                I_new_save(I_dif > 0) = ceil(I_new_save(I_dif > 0));
                I_new_save            = uint8(I_new_save);
                imwrite(I_new_save,save_path,'png');
                I_new                 = double(imread(save_path));
         end
        %comment the following code to speed up 
%         [~,~,keypoints_removed] =  vl_sift(single(I_new),'Octaves',par.octaves ,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',par.cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
%         key_removed(loo, 2)     =  size(keypoints_removed,2);
%         fprintf('removed %d keypoint in this %d th loop', key_removed(loo, 1) - key_removed(loo, 2) , loo);
         psnr                    =  psnrfun(I, I_new, [0, 0])
         key_removed(loo, 3)     =  psnr;
%         max_error_to_I          =  max(max(abs(I-I_new)));
%         key_removed(loo, 4)     =  max_error_to_I;

        %%%%%%%%%%%%%%%%matching%%%%%%%%%%%%%%%%%%%%
        [f_new,d_new]      =   vl_sift(single(I_new),'Octaves',par.octaves ,'Levels',3, 'FirstOctave',par.o_min,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold) ;
        [matches, scores]  =   vl_ubcmatch(d_old, d_new) ;
        [r_Matches]        =   JudgeMat(f_old, f_new, matches);
        key_removed(loo, 5)=   (size(f_old, 2)-size(r_Matches, 2)) / size(f_old, 2);
        fprintf('\ncurrent loo=%d.%d\n, remove rate=%.3f\n',loo,o_loop,key_removed(loo, 5));
        key_removed(loo, 7)=   size(r_Matches, 2);
        %%%%%%%%%%%%%%%end matching%%%%%%%%%%%%%%%%%%%%
        I_cur    =   I_new;
   end                                                      %end octave loop
   loo           =   loo +  1;
end
end


function [keypoints]  =  refine_remove_keys(keypoints, injec_locs)
    ncount  = size(keypoints,1);
    ex_indx = [];
    dis     = 1  ;
    for i = 1: ncount 
        cur_keypoints = keypoints(i,:);
        difx = abs(injec_locs(:,1) - cur_keypoints(1));
        dify = abs(injec_locs(:,2) - cur_keypoints(2));
        ind = difx<=dis & dify<=dis;
        if sum(ind)>0
            %ex_indx= [ex_indx,i];
            fprintf('this should be impossible\n the injected location should be far away from righted matched\n');
        end
    end
    keypoints(ex_indx,:) = [];
end
