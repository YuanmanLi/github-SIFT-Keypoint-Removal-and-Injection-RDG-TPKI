function  f = obj_fun_inject(x, I_col,I,  G_cell, R,win, win_var,var_indx_r, var_indx_c,key_Max, par)
%f = norm(R_current*x-R_current*I);
f = norm(x-I_col);
