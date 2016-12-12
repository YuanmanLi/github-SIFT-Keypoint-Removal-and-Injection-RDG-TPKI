%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Created by Li Yuanman
%%yuanmanx.li@gmail.com
%% This is the implementation of the paper [submitted]: 
%%"SIFT Keypoint Removal via Directed Graph Construction for Color Images" 
%% The code also provides the component of injection algorithm TPKI. 
%% Use the code, please cite the following paper£º
%%Y. Li, J. Zhou, A. Cheng, X. Liu, and Y. Y. Tang, 
%%¡°SIFT keypoint removal and injection via convex relaxation,¡± IEEE TIFS, vol. 11, no. 8, Aug 2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 function [] = main()
   clear all;clc;     % dbstop if all error
   addpath('anti\detector', 'anti\remove','anti\inject'); 
   mkdir('save','inj_save'); mkdir('save','re_save')
   addpath(genpath('vlfeat-0.9.18-bin'));
   current_DB = 'dataset1';
    %%%%%%%%%%%%%%%%%%intitialize parallel  workers%%%%%%%%%%%%%%
    %%%%%%for running over a database%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%     if matlabpool('size')>0
%         matlabpool close;
%     end
%     matlabpool open local 1;    
%     num_worker  =     matlabpool('size');
    num_worker  =     1; 
    in_loop     =     1;
    begin_i     =     6;
    end_i       =     6;         %6 Lena
    N           =     end_i - begin_i+1;
    out_loop    =     ceil(N/ (num_worker*in_loop));
    %%%%%%%%%%%%%%%%%%end intitialization%%%%%%%%%%%%%%%%%%%%  
    
   %%%%%%paras%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   mul_o          =   0;         %do multiple octaves?
   bool_injection =   1;         %do injection? 
   
for  i = 1 : out_loop
    j_min = begin_i + (i-1)*num_worker*in_loop;
    j_max = begin_i + min(i * num_worker*in_loop,N)-1;   
   for im_idx =j_min:j_max       %parfor if parallel
        filename          =  strcat(current_DB,'\',num2str(im_idx),'.pgm');
        dataDB_filename   =  strcat(current_DB,'_',num2str(im_idx));
        mkdir('save\re_save',dataDB_filename);
        re_savepath       =  strcat('save\re_save\', dataDB_filename,'\');
        inj_savepath      =  strcat('save\inj_save\', dataDB_filename,'\');
        I                 =  double( imread(filename));
        utilPara          =  set_Para(mul_o);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        injec_loc = [0 0];  
        %%%%%%%%%%%%%%%%%do remove%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('begin to remove...\n');
        I_cur                 =    I;
        [I_re, key_removed]   =    main_remove(I,I_cur, injec_loc,utilPara,re_savepath,im_idx);
        parsave(strcat('save\',dataDB_filename,'.mat'), key_removed)
        
        if bool_injection == 1
            %%%%%%%%%%%%%%%%do injection%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('begin to inject...\n');
            mkdir('save\inj_save',dataDB_filename);
            [I_inj, injec_loc]  =  main_injection(I,I_re,injec_loc,utilPara,inj_savepath,im_idx);

            %%%%%%%%%%%%%%%%do fine-remove%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf('begin to fine-remove...\n');
            utilPara.save_im    =    0; 
            utilPara.max_loo    =   utilPara.max_loo-30; %less iteration
            [I_refine, key_removed_refine]  =  main_remove(I,I_inj, injec_loc, utilPara,inj_savepath,im_idx);
            parsave(strcat('save\',dataDB_filename,'_refine') ,'key_removed_refine')
            psnr   =  psnrfun(I, I_refine, [0, 0]);
            imwrite(uint8(I_refine),['save\',dataDB_filename, '_re_inj.png']);
        end
         imwrite(uint8(I_re),['save\',dataDB_filename, '_re.png']);
    end
 end
 end
 
function utilPara =  set_Para(mul_o)
    utilPara.peak_threshold        =   4;       
    utilPara.edge_threshold        =   10;
    utilPara.zero_threshold        =   0.001;       
    if mul_o == 1
        utilPara.octaves           =   4;        
        utilPara.o_min             =   -1;
        utilPara.re_octaves        =   [0 3];     % -1 octave is indexed as 0....  
        utilPara.max_loo           =   50;        %set it larger than 80 to get higher KRR
    elseif mul_o==0
        utilPara.octaves           =   1;        
        utilPara.o_min             =   0;
        utilPara.re_octaves        =   [1 1];     % -1 octave is indexed as 0.... 
        utilPara.max_loo           =   40; 
    end
    utilPara.boolCM                =    0;        %do copy move attacking?
    utilPara.save_im               =    1;   
    %%%%%%%%%%%%%%%para for injection%%%%%%%%%%%%%%%%%%%%%%
    if mul_o == 1
        utilPara.inj_octaves       =   [0 1];       % -1 octave is indexed as 0....     
    elseif mul_o ==0
        utilPara.inj_octaves       =   [1 1];     
    end 
    utilPara.minDog_inject         =    3;        %only dog value bigger than this threshold are chosen as candidate location
    utilPara.rank                  =    6;
    utilPara.max_dis_corner        =    2;        %there should at least be one coner in max_dis*max_dis square centered by loc_injection
end

function parsave(fname,data)
    var_name=genvarname(inputname(2));
    eval([var_name '=data;']);
     save(fname,var_name);
end