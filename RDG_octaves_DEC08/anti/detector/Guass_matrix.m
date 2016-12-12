%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function P = Guass_matrix(sigDiff, I)
[x]  =  2*max(ceil(4.0*sigDiff),1)+1;
Gauss_matrix = fspecial('gaussian', [x x], sigDiff);  

win  = size(Gauss_matrix);
[w h]=size(I);
bounder=[floor(win(1)/2),floor(win(2)/2)] ;

idx_pic_m = reshape(1:w*h, w,h);
idx_pic_m = padarray(idx_pic_m, bounder, 'replicate'); 
%%we only calculated a patch located in a block. only patch are variable,
%%others are constrant. 

r = bounder(1)+1:w+bounder(1);
c = bounder(2)+1:h+bounder(2);
n = length(r)*length(c);
non_zeros_indx = zeros(n, win(1)*win(2));

k = 0;
for i =1: win(1)
    for j =1: win(2)
            k = k+1;
            temp_idx =  idx_pic_m(r-bounder(1)+i-1, c-bounder(2)+j-1);
            temp_idx = temp_idx(:);
            non_zeros_indx(:,k) = temp_idx;
    end
end

row_vector = repmat(1:n, [win(1)*win(2), 1]);  % generate the row index with nonzero elements
row_vector = row_vector(:);
non_zeros_indx = non_zeros_indx';          % generate the column index with nonzero elements
non_zeros_indx = non_zeros_indx(:);
temp_G = Gauss_matrix';
val    = repmat(temp_G(:), n,1);                  % generate the nonzero values corresponding to row and column
P      = sparse(row_vector, non_zeros_indx, val,n,n, length(non_zeros_indx));
 