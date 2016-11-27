% compile_filter

disp('compiling/installing adaptivefilt_mex...')
mex adaptivefilt_mex.c
eval(['!move adaptivefilt_mex.' mexext ' ../']);

