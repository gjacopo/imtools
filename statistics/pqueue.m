%% PQUEUE - Implement a priority queues.
%
%% Syntax
%     h = PQUEUE('init');
%     h = PQUEUE('insert', hin, v);
%     [h,v] = PQUEUE('remove', hin );
% 
%% Inputs
% *|action|* : string specifying the action to perform on a priority queue,
%      it is either |'init'|, |'insert'| or |'remove'|.
%
% *|hin|* : structure representing a priority queue (aka, a min/max heaps),
%      with fields |tree| and |count|, where to insert or remove an element.
%      
% *|v|* : value to insert into the queue |hin| when |action='insert'|.     
% 
%% Outputs
% *|h|* : updated priority queue.
%
% *|v|* : when |action='remove'|, variable removed. 

%% Function implementation
% -------------------------------------------------------------------------
function [h, varargout] = pqueue(action, varargin)

if (strcmp(action,'init') && nargin~=0) || ...
        (strcmp(action,'remove') && nargin~=1) || ...
        (strcmp(action,'insert') && nargin~=2)
    error('pqueue:inputerror','unmatched inputs with selected action');
end

if nargin>=2,  hin = varargin{1};
    if nargin==3,  val = varargin{2};  end
end

%%
% implements action on priority queue 

switch action
    case 'init'
        h = heapinit();
    
    case 'insert'
        h = heapinsert(hin,val);
        
    case 'remove'
        [val,h] = heapremove(hin);
        if nargout==2,  varargout{1} = val;  end
end

end % end of pqueue


%% Subfunctions

%%
% |HEAPINIT| - Return a heap h that is empty.
% This must be called before the heapinsert and heapremove routines.
% -------------------------------------------------------------------------
function h = heapinit()
h.count = 0;
h.tree = [];
end


%%
% |HEAPINSERT| - Return the heap created by inserting the value val into the
% existing heap hin.
% -------------------------------------------------------------------------
function hout = heapinsert(hin,val)

if (hin.count == 0)
    hout.count = 1;
    hout.tree = val;
else
    hout.count = hin.count+1;
    hout.tree = [hin.tree val];
end

cur = hout.count;
parent = floor(cur/2);
found = 0;

while (found == 0)
    if (parent == 0)
        found = 1;
    else
        if (hout.tree(parent) > hout.tree(cur))
            tmp = hout.tree(parent);
            hout.tree(parent) = hout.tree(cur);
            hout.tree(cur) = tmp;
            cur = parent;
        else
            found = 1;
        end
    end
    
    parent = floor(cur/2);
end
end


%%
% |HEAPREMOVE| - Remove the root element of the heap, and return it and the
% new heap with the element removed.
% -------------------------------------------------------------------------
function [v,h] = heapremove(hin)

h = hin;
v = h.tree(1);
h.tree(1) = h.tree(h.count);
h.count = h.count - 1;
if (h.count == 0)
    h.tree = [];
else
    h.tree = h.tree(1:h.count);
end

cur = 1;
lchild = 2;
rchild = 3;
found = 0;

while (found == 0)
    numchildren = (lchild <= h.count) + (rchild <= h.count);
    
    if (numchildren == 0)
        found = 1;
    elseif (numchildren == 1)
        if (h.tree(lchild)<h.tree(cur))
            tmp = h.tree(lchild);
            h.tree(lchild) = h.tree(cur);
            h.tree(cur) = tmp;
            cur = lchild;
        else
            found = 1;
        end
    else
        [tmp, idx] = min([h.tree(lchild) h.tree(rchild)]);
        if (idx == 1)
            idx = lchild;
        else
            idx = rchild;
        end
        
        if (tmp < h.tree(cur))
            h.tree(idx) = h.tree(cur);
            h.tree(cur) = tmp;
            cur = idx;
        else
            found = 1;
        end
    end
    
    lchild = cur*2;
    rchild = lchild+1;
end

end