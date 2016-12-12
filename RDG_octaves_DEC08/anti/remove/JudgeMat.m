%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% used to judge whether this match pair is wrong or not
function [r_Matches] = JudgeMat(f_old, f_new, matches)
    max_dis  = 8;
    key_old  = f_old(1:2,:);
    key_new  = f_new(1:2,:);
    matched_old = key_old(:,matches(1,:));
    matched_new = key_new(:,matches(2,:));
    difx = abs(matched_old(1,:) - matched_new(1,:) );
    dify = abs(matched_old(2,:) - matched_new(2,:) );
    dis = sqrt(difx.^2+dify.^2);
    r_Matches_idx = dis<max_dis;
    r_Matches     =  matches(:,r_Matches_idx);
    [~,r1]        =  size(matches); [~,r2] = size(r_Matches);
    fprintf('before fine-matching %d, after fine-matching %d\n', r1, r2);
end
% r_Matches = [];
% for i = 1 : size(matches, 2)
%     x_old = f_old(1, matches(1,i)) ;        
%     x_new = f_new(1, matches(2,i)) ; 
%     y_old = f_old(2, matches(1,i)) ;  
%     y_new = f_new(2, matches(2,i)) ;
%     r_x1 = x_old - 5; % error range is 5
%     r_y1 = y_old - 5;
%     r_x2 = x_old + 5;
%     r_y2 = y_old + 5;
%     new = [x_new, y_new];
%     r1 = [r_x1, r_y1];
%     r2 = [r_x2, r_y2];
%     dif1 = new - r1;
%     dif2 = new - r2;
%     if (sum(dif1(1, :) > 0) == 2 && sum(dif2(1, :) < 0) == 2)
%         r_Matches = [r_Matches, matches(:, i)];
%     end
% end

