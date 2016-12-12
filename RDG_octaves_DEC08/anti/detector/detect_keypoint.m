%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [keypoints,G_layers] = detect_keypoint(I, par)
win      = [3,3];                                                 %
zero_thr = par.zero_threshold;
[h w]    = size(I);
keypoints= [];
%%%%%%%%%%parameter setting %%%%%%%%%%%%
f.nlevels = 3;
f.octave  = 1;
f.o_min   = 0;
f.s_min   = -1;

f.sigman = 0.5;
f.sigmak = 2.0^(1.0/f.nlevels);
f.sigma0 = 1.6 * f.sigmak;
f.dsigma0 = f.sigma0 * sqrt(1.0 - 1.0/(f.sigmak*f.sigmak));

f.peak_thresh = par.peak_threshold ;                                        %this is used for finding local minma or maxma
f.edge_thresh = par.edge_threshold ;
f.norm_thresh = 0;
f.magnif      = 3;
f.windowSize  = 2;

if f.o_min<0
    error 'handle upsample'
elseif f.o_min>0
    error 'handle downsample'
else
    I_resize = I;                                                 %just the original image
end
[a b]  = size(I_resize);
detect_size = a * b;                         %every loop, just detector a subimage of size
sa     = f.sigma0* (f.sigmak^ f.s_min);
sb     = f.sigman * 2.0^(-f.o_min);
loops  = ceil(w*h/detect_size);
cols_loop = ceil(detect_size/h);             %cols to detect in every loop
%%handle the first level
if(sa > sb)
    sd = sqrt(sa*sa - sb*sb);
end
%guassian.sigma = sd;
G_cell       = cell(1,6);
Sig_cell     = zeros(1,6);
Sig_cell(1)  = sd;
[x]          =getGuassianSize(sd);
padding(1)   = floor(x/2)+2;
idx = 1;
for  s = f.s_min+1 : f.nlevels - f.s_min        %calculate all the sigmas
    idx = idx+1;
    sd = f.dsigma0 * (f.sigmak^ s) ;
    [x]=getGuassianSize(sd);
    padding(idx) = floor(x/2)+2;
end
pad = max(padding);
fprintf('image is divided into %d segments\n', loops);
for Loo = 1:loops
    fprintf('Loo=%d  ', Loo);
    if(Loo ==1)
        I_resize = I(:, 1:min( cols_loop + pad,w));
    elseif(L==loops)
        I_resize = I(:, max(1,(Loo-1)*cols_loop-pad+1): end);
    else
        I_resize = I(:, max((Loo-1)*cols_loop-pad+1,1): min(Loo*cols_loop+pad, w));
    end   
    [x]= getGuassianSize(Sig_cell(1));
    G_cur = fspecial('gaussian', [x x], Sig_cell(1));
    G_cell{1} = G_cur;
    idx = 1;
    for  s = f.s_min+1 : f.nlevels - f.s_min
        idx = idx+1;
        sd  = f.dsigma0 * (f.sigmak^ s) ;
        Sig_cell(idx) = sd;
        %G_cur =  Guass_matrix(sd, I_resize);
        [x] = getGuassianSize(Sig_cell(idx));
        G_cur = fspecial('gaussian', [x x], Sig_cell(idx));
        G_cell{idx} = G_cur;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%generate the Guassian matrice%%%%%%%%%%%%%%
    G1_I = imfilter(I_resize,G_cell{1},'replicate');
    G2_I =  imfilter(G1_I,G_cell{2},'replicate');
    G3_I =  imfilter(G2_I,G_cell{3},'replicate');
    G4_I =  imfilter(G3_I,G_cell{4},'replicate');
    G5_I =  imfilter(G4_I,G_cell{5},'replicate');
    G6_I =  imfilter(G5_I,G_cell{6},'replicate');
    G_layers = {G1_I,G2_I,G3_I,G4_I,G5_I,G6_I};
    
    D1 =  G2_I-G1_I;
    D2 =  G3_I-G2_I;
    D3 =  G4_I-G3_I;
    D4 =  G5_I-G4_I;
    D5 =  G6_I-G5_I;
    
    Dog = {D1, D2, D3, D4, D5};
    
    r   = 2:a-1;                                               %the boudary can't be compared
    c   = 2:b-1;
    keyPointer_layer = cell(1,3);
    for L= 1:3                                                    %three layers
        patches = zeros(27, (a-2)*(b-2));
        k = 1;
        for z = L:L+2
            for i = 1:win(1)
                for j = 1:win(2)
                    Dz   = Dog{z};
                    temp = Dz(r+i-2, c+j-2);
                    temp = reshape(temp',  1, (a-2)*(b-2));   %here has a transposition. then ordered as row
                    patches(k, :) = temp;
                    k = k + 1;
                end
            end
        end
        center  = repmat(patches(14,:), 27,1);
        patches_dif = center - patches;                % if is keypoint, then also be less than zero, or bigger than zero
        patches_dif(abs(patches_dif) <= zero_thr) = 0;
        patches_dif_sign  = sign(patches_dif);
        patches_dif_sum   = sum(patches_dif_sign);
        key_idx           = find(abs(patches_dif_sum)==26); %get the index of the keypoints in patches
        
        %%vl_feat refines the constrast keypoints
        cent_Dog            =  patches(14,:);
        cand_key            =  cent_Dog(key_idx);
        constrast_refined   =  key_idx(cand_key > 0.8 * f.peak_thresh | cand_key< -0.8 * f.peak_thresh);
        %%end refine%%%%%%%%%%%%%%%%%%%%%%
        
        rows                = floor((constrast_refined-1) / (b-2))+2;
        cols                = mod((constrast_refined-1), b-2) +2;
        keyPointer_layer{L} = [rows', cols', L*ones(length(rows),1)];
    end
    keypoints_Loo = union(keyPointer_layer{1}, keyPointer_layer{2}, 'rows');
    keypoints_Loo = union(keyPointer_layer{3}, keypoints_Loo, 'rows');
    
    if(Loo ==1)
        idx = find(keypoints_Loo(:,2,:) > cols_loop);
        sub_val = 0;
    elseif(Loo==loops)
        idx = find(keypoints_Loo(:,2,:) < pad+1);
        sub_val = pad;
    else
        idx = find(keypoints_Loo(:,2,:) < pad+1 | keypoints_Loo(:,2,:) > cols_loop+pad);
        sub_val = pad;
    end
    keypoints_Loo(idx,:,:) = [];
    keypoints_Loo(:,2,:)   = keypoints_Loo(:,2,:)+(Loo-1)*cols_loop-sub_val;
    keypoints = [keypoints; keypoints_Loo];
end
%save 'keypoints.mat'  keypoints ;



