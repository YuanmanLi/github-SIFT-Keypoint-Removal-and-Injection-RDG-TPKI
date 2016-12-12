%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [refined_loc] = refine_candidates(candidate_Loc,key_loc,coner_loc,para)
    
    max_dis  =  para.max_dis;   %%%%there should at least be one coner in max_dis*max_dis square centered by loc_injection
    min_dis  =  para.min_dis;    %%%%there are not any orignal keypoints in min_dis*min_dis square centered by loc_injection
    ncount   =  size(candidate_Loc,1);
    refined_loc = [];
    
    %%%%%%%%%%map the keypoint location to orginal image %%%%%%%%%%%%%%
    temp1 = candidate_Loc(:,1);   temp2 = candidate_Loc(:,2); 
    temp1 = ceil(temp1*(2^(para.cur_o-1))); temp2 = ceil(temp2*(2^(para.cur_o-1)));  
    remap_Loc = [temp1,temp2];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i     =  1:ncount
        cur_cand      =   remap_Loc(i,1:2);
        dif_cx        =  abs(coner_loc(:,1) - cur_cand(1));
        dif_cy        =  abs(coner_loc(:,2) - cur_cand(2));
        indx          =   dif_cx<=max_dis & dif_cy <=max_dis;
        if(sum(indx)<=0)
            continue;                                 %current location is not near a coner
        end
        dif_keyx      =  abs(key_loc(:,1) - cur_cand(1));
        dif_keyy      =   abs(key_loc(:,2) - cur_cand(2));
        dis = sqrt(dif_keyx.^2+dif_keyy.^2);                   %this distance method is from TIFS
        indx          =   dis<=min_dis;
        if (sum(indx)>0)          
            continue;       %current location is  near some original keypoints
        end
        refined_loc = [refined_loc; candidate_Loc(i,1:2)];
    end

end