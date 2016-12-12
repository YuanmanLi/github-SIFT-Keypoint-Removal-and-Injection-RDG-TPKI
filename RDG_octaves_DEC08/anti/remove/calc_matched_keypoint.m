%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matched_keypoints] =  calc_matched_keypoint(keyrep, r_match,octaves)
%[~,octaves]       = size(keyrep);
nkeys     = keyrep(end-octaves+1:end);
keyrep      = keyrep(1:end-octaves);
start_loc   = 1; 


matched_keypoints = cell(octaves,1);
r_match_new       = r_match(2,:);                        %get the index of righted match decriptor
nkeys_o           = zeros(octaves,1);
cum_descr = 0;
for i = 1: octaves
    cur_keys    = keyrep(start_loc: start_loc+nkeys(i)-1);
    start_loc   = start_loc + nkeys(i);

    nkeys_o(i) = nnz(cur_keys);
    n_descr    = sum(cur_keys);                      %number of descrptors in the current octave
    map_key_descr = zeros(n_descr,1);
    
    
    %%%%%%%%%%%%%%%map keypoint and descriptor%%%%%%%%%%%%
    j = 1; 
    for  k= 1: nkeys_o(i)                                  
        temp = cur_keys(k);
        map_key_descr(j) = k;
        while temp-1 > 0
             map_key_descr(j+1) =map_key_descr(j) ;
             j    = j+1;
             temp = temp-1;
        end
        j = j+1;
    end
    if j ~= n_descr+1
        error('map_key_descr error!!\n');
    end
    temp_idx     = (cum_descr+n_descr>=r_match_new & r_match_new>cum_descr);
    r_descr_indx = r_match_new(temp_idx) - cum_descr;           %find right match indx of descriptor                     
    r_key_indx   = map_key_descr(r_descr_indx);                         %find right match indx of keypoints 
    matched_keypoints{i} =  r_key_indx;
    cum_descr    = cum_descr + n_descr;
end
end