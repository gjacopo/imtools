% compile_pyramid

disp('compiling/installing reduced2d_mex...')
mex -I. messagedisplay.c pyramidfilters.c pyramidtools.c reduce2d_mex.c ...
    -output reduce2d_mex
eval(['!move reduce2d_mex.' mexext ' ../']);


disp('compiling/installing expand2d_mex...')
mex -I. messagedisplay.c pyramidfilters.c pyramidtools.c expand2d_mex.c ...
    -output expand2d_mex
eval(['!move expand2d_mex.' mexext ' ../']);

