% compile_feature.m

disp('compiling/installing anisoedge_mex...')
mex anisoedge_mex.c
eval(['!move anisoedge_mex.' mexext ' ../']);

disp('compiling/installing rothwelledge_mex...')
mex rothwelledge_mex.c
eval(['!move rothwelledge_mex.' mexext ' ../']);

disp('compiling/installing susan_mex...')
mex susan_mex.c
eval(['!move susan_mex.' mexext ' ../']);

disp('compiling/installing fast_mex...')
mex -I./fast-C-src-2.1/ fast-C-src-2.1/fast_11.c fast-C-src-2.1/fast_9.c ...
    fast-C-src-2.1/fast_10.c fast-C-src-2.1/fast_12.c fast-C-src-2.1/nonmax.c ...
    fast-C-src-2.1/fast.c fast_mex.c -output fast_mex
eval(['!move fast_mex.' mexext ' ../']);
