%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%% This is the implementation of the paper: "SIFT Keypoint Removal via Directed Graph Construction for Color Images"
%% The code also provides the component of injection algorithm TPKI. 
%% Use the code, please cite the following paper£º
%% Y. Li, J. Zhou, A. Cheng, X. Liu, and Y. Y. Tang, 
%% ¡°SIFT keypoint removal and injection via convex relaxation,¡± IEEE TIFS, vol. 11, no. 8, Aug 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [I_cur, injected_key]     =     main_injection(I,I_cur,injected_key,utilPara,inj_savepath,im_idx)
%%%%%%paras%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.zero_threshold       =    utilPara.zero_threshold ;    
par.peak_threshold       =    utilPara.peak_threshold;                                      %parameter for vl_feat constrast¡£ used in detect_keypoint.m
par.edge_threshold       =     utilPara.edge_threshold;  
par.minDog_inject        =     utilPara.minDog_inject;    %only dog value bigger than this threshold are chosen as candidate location

par.rank       =  utilPara.rank;                                  %only location with dog in the topper rank are considered to inject
par.o_min      = utilPara.o_min;
par.octaves    = utilPara.octaves ;     
par.inj_octaves= utilPara.inj_octaves;


par.con_threthold   = 0.1;  
par.cur_s          =  1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%getting the keypoints of all the ocataves%%%%%%%%%%%%%%%%%%%%%%%%%%
key_loc = get_all_octaves_keypoints(I, par); %keypoint location of all the octaves 

%%%%%%%%%%%%%%%%%%%%find the corner%%%%%%%%%%%%%%%
%%%Shi & Tomasi's minimum eigenvalue method%%%%%%%%%%%%%%%%
N_strongest_corner = 5000;
Corners             =   detectMinEigenFeatures(I,'FilterSize',3,'MinQuality',0.00001);
strCornes          =    Corners.selectStrongest(N_strongest_corner);
coner_loc          =    round(double(strCornes.Location));
%%%%%%%%%%%%get candidate location in specific octave%%%%%%%%%%%
bool_break = 0;
iter = 1;
while iter < 10
  for i = par.inj_octaves(1):par.inj_octaves(2)
      for j = 1:3
          par.cur_o       = i;
          par.cur_s = j;
          if par.cur_o == 0
             par.con_threthod1        = 0.1;   
             par.varwin               = [3  3]; 
             par.stop_max_error       = 100;  
          elseif par.cur_o ==1
             par.con_threthod1        = 0.05;   
             par.varwin               = [5  5]; 
             par.stop_max_error       = 150;  
          else
             par.con_threthod1        = 0.03;   
             par.varwin               = [5  5]; 
             par.stop_max_error       = 150;  
          end
         par.stop_max_error =  par.stop_max_error+ 5*(iter-1);
         [~,candidate_Loc]   =  Get_Candidates_Loc(I_cur, par);     
         %%%%%%%%%%%%refine the cadidate injection locations%%%%%%%%%%% 
         refinePara. max_dis   = utilPara.max_dis_corner;   %%%%there should at least be one coner in max_dis*max_dis square centered by loc_injection                                                              %%%because original not every                                                                    %%%fine-remove process
         refinePara. min_dis  = 8;    %%%%there are not any orignal keypoints in min_dis*min_dis square centered by loc_injection
         refinePara.cur_o     = par.cur_o;
         refined_loc          = refine_candidates(candidate_Loc{par.cur_s},key_loc(:,1:2),coner_loc,refinePara);
%          tic
         [I_new,injected_key] = inject_keypoint(I,I_cur,refined_loc,injected_key, par);
%          remove_time = toc; 
%          fprintf('in this interation, time = %s\n', remove_time); 
         I_cur = I_new;  
          
          %%%%added 25-Aug-2015%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          [f_new,d_new,keypoints,keyrep]    =           vl_sift(single(I_new),'Octaves',par.octaves,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',par.cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
          %psnr                                                      =           psnrfun(I, I_new, [0, 0])
          [f_new1,d_new1,keypoints1,keyrep] =           vl_sift(single(I),'Octaves',par.octaves,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',par.cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
          [matches, scores]                 =           vl_ubcmatch(d_new1,d_new) ;
          [r_Matches]                       =           JudgeMat(f_new1, f_new, matches);
          save_path                         =            strcat(inj_savepath,int2str(im_idx),'_inj_iter', int2str(iter), '.png');
          imwrite(uint8(I_new),save_path);
          if size(d_new,2) - size(d_new1,2)-0.8*size(r_Matches,2) >0 
              bool_break =1;   %inject enough
              break;       
          end  
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      end
      %%%%after each octave, exclude injected location at which the
      %%%%injected keypoint dissmissed%%%%%%%%%%%%
      %injected_key = refine_injec_locs(injected_key, I_new,par);
       %%%%added 25-Aug-2015%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if bool_break == 1
        break
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end
  
  save_path    =   strcat(inj_savepath,int2str(im_idx),'_inj_iter', int2str(iter), '.png');
  imwrite(uint8(I_new),save_path);
  if bool_break ==1
      break;  
  else
     iter = iter +  1;
     par.minDog_inject   = par.minDog_inject   - 0.4;
     N_strongest_corner = N_strongest_corner + 1000;
     Corners             =   detectMinEigenFeatures(I,'FilterSize',3,'MinQuality',0.00001);
     strCornes          =    Corners.selectStrongest(N_strongest_corner);
     coner_loc          =    round(double(strCornes.Location));
  end
end