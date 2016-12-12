%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function key_loc = get_all_octaves_keypoints(I_cur, par)
key_loc = [];
for i= par.o_min: par.o_min + par.octaves-1
    par.cur_o                                         =           i + 1;
    [~,~,keypoints,~]       =           vl_sift(single(I_cur),'Octaves',par.octaves,'Levels',3, 'FirstOctave',par.o_min,'outputOctave',par.cur_o-1,'ThreEqual',par.zero_threshold,'PeakThresh', par.peak_threshold, 'EdgeThresh', par.edge_threshold);
    keypoints               =           keypoints';  
    [~,idx]                 =           sort(keypoints(:,1));      
    keypoints               =           keypoints(idx,:);       
    
    [keypoints1]            =           detect_keypoint_nextO(I_cur, par);
    keypoints               =           fit_mat_C(keypoints, keypoints1);
    key_loc                 =            [key_loc; keypoints];       
end