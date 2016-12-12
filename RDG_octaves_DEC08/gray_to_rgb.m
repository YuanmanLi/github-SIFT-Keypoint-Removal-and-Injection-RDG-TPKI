%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%do not use w = [0.2989, 0.5870, 0.1140] directly. lost precision%%%%%
clc 
clear all
dbstop if error

ori_color_im = '';
re_gray_im = ''; 
save_path = [re_gray_im(1:end-4),'_color_recon.png'];

T = inv([1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703]);
w = T(1,:);


im_rgb = double(imread(ori_color_im));
im_rgb_r = im_rgb(:,:,1);  im_rgb_g = im_rgb(:,:,2); im_rgb_b = im_rgb(:,:,3); 

im_gray = double(rgb2gray(uint8(im_rgb)));
im_gray_hat =double(imread(re_gray_im));

im_dif = double(im_gray) - im_gray_hat;
x_indx = find(im_dif ~= 0);
Y1 = [im_rgb_r( x_indx), im_rgb_g( x_indx), im_rgb_b( x_indx)]';
%%%%%%%%%%%for test%%%%%%%%%%%%%%%%%%%%
im_gray_new = double(uint8(w*Y1));
aa = im_gray_new' - im_gray(x_indx);  %%aa should all be zero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

invAAt = 1/(w*w');

%%%%%%%%%%%%%%%%%strategy 1, slow%%%%%%%%%%%%%%%%%
% Y = double(Y1(:));
% mask = repmat({w}, length(x_indx),1);
% A = blkdiag(mask{:});
% C = double(im_gray_hat(x_indx));
% AY_C = A*Y-C;
% X_temp = Y - invAAt* A' * AY_C;
% X = uint8(reshape(X_temp, 3,  length(x_indx)));
% 
% %%%%%%%%%%%%%%%%stratege 2 fast%%%%%%%%%%%%%%%%%%%%%%
C = double(im_gray_hat(x_indx));
X = Y1 - invAAt * w'*(w*Y1 - C');

%%%%%%%%%some values are bigger than 255 or smaller than 0%%%%%%%%%%%%%
X = uint8(X);
dif_comp = round(w*double(X) - C');
x_indx_comp = find(dif_comp~=0);
N_comp = length(x_indx_comp);

max_vars = 100;
if N_comp > max_vars
    fprintf('the values of more than 1000 pixels are near 0 or 255. suggest split variables into small block to make lsqlin run fast\n');
end
loops = ceil(N_comp/max_vars);
x_indx_comp_pre =  x_indx_comp; 
for i = 1 : loops 
        fprintf('loop %d/%d  ', i, loops);
        indx= (i-1)*max_vars+1: min(i*max_vars, N_comp);
        x_indx_comp =  x_indx_comp_pre(indx);
        n_vars = length(indx);

        Y_comp =Y1(:,x_indx_comp); 
        Y_comp = double(Y_comp(:));
        mask = repmat({w}, length(x_indx_comp),1);
        A_eq = blkdiag(mask{:});
        C_eq = C(x_indx_comp);
        X_comp = lsqlin(eye(3*n_vars), Y_comp,[],[], A_eq, C_eq, 0*ones(1,3*n_vars), 255*ones(1,3*n_vars));
        X_comp = reshape(X_comp, 3, n_vars);

        %%%%%%%%maximum error would be 1 because the round error%%%%%%%%
        dif_comp_2 = round(w*double(uint8(X_comp))-C_eq');
        x_indx_comp_2 = find(dif_comp_2~=0);

        if ~isempty(x_indx_comp_2)
            X_comp  = uint8(repmat(-1*sign(dif_comp_2), 3,1).*ones(3, n_vars) + double(X_comp));
            dif_comp_2 = round(w*double(uint8(X_comp))-C_eq');
            x_indx_comp_2 = find(dif_comp_2~=0);

            if ~isempty(x_indx_comp_2)       %here should exactly recover the gray-scale image
                error('this can not happen!\n');
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X(:,x_indx_comp) = uint8(X_comp);
end


im_recon_r = im_rgb_r; im_recon_g = im_rgb_g; im_recon_b = im_rgb_b;
im_recon_r(x_indx) = X(1,:); 
im_recon_g(x_indx) = X(2,:); 
im_recon_b(x_indx) = X(3,:); 
im_recon(:,:,1) = im_recon_r; im_recon(:,:,2) = im_recon_g; im_recon(:,:,3) = im_recon_b; 
im_recon =  uint8(im_recon);
imshow(im_recon);
imwrite(im_recon, save_path);




%%%%%%%%%%%validate%%%%%%%%%%%%%%%%%%
im_gray_recon = double(rgb2gray(im_recon));
im_dif_recon = im_gray_recon - im_gray_hat;
