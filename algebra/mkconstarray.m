function out = mkconstarray(class, value, size)
%MKCONSTARRAY creates a constant array of a specified numeric class.
%   A = MKCONSTARRAY(CLASS, VALUE, SIZE) creates a constant array 
%   of value VALUE and of size SIZE.

%   Copyright 1993-2000 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2000/01/21 20:19:08 $

out = repmat(feval(class, value), size);

