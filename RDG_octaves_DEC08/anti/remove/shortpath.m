%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Created by Li Yuanman
%% yuanmanx.li@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [path, D ,U, U_pre, b_flag]  = shortpath(G, start_indx, key_Max)

% G =  [ inf  6   1      12     inf;
%             -6  inf  inf   2       inf;
%             -1  inf     inf   inf    9;
%          -12  -2     inf      inf    -1;
%           inf   inf    -9      1       inf]
% 
%  start_indx = 1;
%  key_Max = -1;
 
% if nargin < 2
%     error('input args are not enough!')
% end
if key_Max == 1
     temp_idx    =  find(G == 10^6);
     G           =  -G;
     G(temp_idx) =  10^6;
end
[a a1]  = size(G);

if (a ~= a1)
    error('the adjacent matrix should be square!')
end
if key_Max ~=-1 && key_Max ~= 1
    error('key_Max should be either -1 and 1')
end
b_flag        =   0;                                           % b_flag ==1 means find a good path which guarantee no pixel being generated
V             =       zeros(1, a);  
V(start_indx) =       1;                                             % the start_indx is marked as selected
temp_p        =       start_indx;                            % the sart_indx is marked as intermediate point
D             =       inf *ones(1, a);                     % the shortest distance from start_indx
D(start_indx) =        0;

U             =       [temp_p];                               %the order of seleting points
U_pre         =        zeros(1,a);                            % the previous point for each point. this is used to inversely find the path
U_pre(1)      =        temp_p;

b_terminal    =  0;
num_p_s       =   a/5;                                                 % 5 means five Dog.  num_p_s means number of pixels(including the boundary) per layer
while  sum(V) ~=   a  && ~b_terminal
    non_sel    =       find(V == 0);
    D(non_sel) =       min (D(non_sel),  D(temp_p) + G(temp_p, non_sel) );           %get the new shortest patch by considering the intermediate point
    new_indx   =       find(D(non_sel) == min(D(non_sel)));                                         % get the new shortest path from start the sart_p to the non_seleted point
    
   
    temp_p     =        non_sel(new_indx(1));                                                                       % the new point to be selected    
    V(temp_p)  =       1;
    U          =       [U temp_p];
    indx       =       find(D(U)' ==  D(temp_p) - G(U, temp_p));
    U_pre(temp_p)    =      U(indx(1));
    
   %%%%%% begin terminal conditions  %%%%%%%%%%%%%%%%%

    if key_Max == -1 && (G(U(indx(1)),temp_p) < 0 | G(U(indx(1)),temp_p) == 10^6 |  G(temp_p,U(indx(1))) == 10^6)
         if G(U(indx(1)),temp_p) < 0 || temp_p <= num_p_s  || temp_p >= num_p_s*4+1     %find the good path, no keypoint will be generated
            fprintf('good short path,')
            b_flag           =  1;
         end
          b_terminal   = 1;
    elseif key_Max == 1 && (G(U(indx(1)),temp_p) < 0 | G(U(indx(1)),temp_p) == 10^6 |  G(temp_p,U(indx(1))) == 10^6)
         if G(U(indx(1)),temp_p) < 0|| temp_p <= num_p_s  || temp_p >= num_p_s*4+1  %find the good path, no keypoint will be generated
             fprintf('good short path,')
            b_flag          =       1;
         end
          b_terminal = 1;
    end
        %%%%%%end   terminal conditions  %%%%%%%%%%%%%%%%%
end

 %%%%%%%%%%%%get the path using U and U_pre 
 end_p      =            temp_p;        
 path       =            end_p;
 temp       =             end_p;
 while temp  ~=  start_indx; 
      temp   = U_pre(U(U == temp));                      % the corresponding point in U_pre is the previous point to end_point
      path   = [temp path];
 end
 fprintf('path:');
 fprintf('%d  ',path);
 fprintf('\n',path);

 