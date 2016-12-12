%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function injec_locs = refine_injec_locs(injec_locs, I_inject,par)
    dis  = 1;                            %because C and matlab have a small margin of error £¨set dis=0 just ignore such phenomenon£©
   %%%%%%%%%%%%%%%get all the keypoints and remap to original size
   keypoints_src = [];
    for i= par.o_min: par.o_min + par.octaves-1
        cur_o                   =           i + 1;
        [~,~,keypoints,~]       =           vl_sift(single(I_inject),'Octaves',par.octaves,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
        keypoints               =           keypoints';  
        [~,idx]                 =           sort(keypoints(:,1));      
        keypoints               =            keypoints(idx,:);       
        %%%%%%%%%%map the keypoint location to orginal image %%%%%%%%%%%%%%
        temp1 = keypoints(:,1);   temp2 = keypoints(:,2); 
        temp1 = ceil(temp1*(2^(cur_o-1))); temp2 = ceil(temp2*(2^(cur_o-1)));  
        remap_Loc = [temp1,temp2];
        keypoints_src = [keypoints_src; remap_Loc];
    end
    
    ncount = size(injec_locs,1);
    ex_indx = [];
    for i = 1: ncount 
        cur_injec_loc = injec_locs(i,:);
        difx = abs(keypoints_src(:,1) - cur_injec_loc(1));
        dify = abs(keypoints_src(:,2) - cur_injec_loc(2));
        ind  = difx<=dis & dify<=dis;
        if sum(ind)<=0
            ex_indx= [ex_indx,i];
        end
    end
    injec_locs(ex_indx,:) = [];