%% FASTCPDA_BASE - Base function for FASTCPDA.
%
%% Syntax
%     [pt, map, cd] = FASTCPDA_BASE(BW, Gap_size, T_angle, EP);
%
%% Acknowledgement
% This is an adaptation of the original code of [AL08,ALFR09] available at:
%  <http://www.mathworks.com/matlabcentral/fileexchange/28207-a-fast-corner-detector-based-on-the-chord-to-point-distance-accumulation-technique>
%
% Part of this code was already from the source code of [HY04,HY08] available at:
%  <http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7652&objectType=file>
%
%% See also
% Related:
% <FASTCPDA.html |FASTCPDA|>,
% <HARRISCORNER_BASE.html |HARRISCORNER_BASE|>,
% <SUSANCORNER_BASE.html |SUSANCORNER_BASE|>,
% <EDGECORNER_BASE.html |EDGECORNER_BASE|>,
% <FASTCORNER_BASE.html |FASTCORNER_BASE|>.

%% Function implementation
function [pt, varargout] = fastcpda_base(BW, Gap_size, T_angle, EP)

[sizex sizey C] = size(BW);

%%
% dealing with multispectral images

pt = cell(C,1);
if nargout>=2 
    varargout{1} = false(size(BW));
    if nargout==3, varargout{2} = cell(C,1); end;
end

if C>1
    for c=1:C
        [tmp1, tmp2, tmp3] = fastcpda_base(BW(:,:,c), Gap_size, T_angle, EP);
        pt{c} = tmp1{1};
        if nargout>=2, 
            varargout{1}(:,:,c) = tmp2;
            if nargout==3, varargout{2}{c} = tmp3{1};  end;
        end
     end
    return;
end

%% 
% extract curves from the edge-image
[curve, curve_start, curve_end, curve_mode, curve_num, TJ, BW_edge] = ...
    extract_curve(BW, Gap_size);  
% curve : MATLAB cell data structure where each cell is a 2D array containing
%         pixel locations (x and y values)
% curve_start : starting points of extracted curves
% curve_end : ending points of extracted curves
% curve_mode : two curve modes - 'line' and 'loop'. If the both ends of
%              a curve are at maximum 25 square pixels (default) away, then 
%              the curve is a loop curve, otherwise a line curve
% curve_num : number of extracted curves
% TJ : the T-junction found in the edge-extraction process
BW_edge;                                                               %#ok

if ~isempty(curve)
    
    [point_selected, smoothed_curve, Width] = ...
        selectPoint(curve, curve_mode, curve_num, sizex, sizey);
    
    % detect corners on the extracted edges
    [pt{1}, ~, cd2] = ...
        getcorner(curve, curve_mode, curve_start, curve_num, T_angle, ...
        point_selected, smoothed_curve, Width);
    % pt : n by 2 matrix containing the positions of the detected corners,
    %      where n is the number of detected corners
    % index : MATLAB cell data structure where each cell is an 1D column
    %         matrix contaning the edge pixel numbers (in curve) where the
    %         corners are detected
    % Sig : the sigma values used to smooth the curves
    % cd2 : cpda curvature values of the detected corners
    
    % update the T-junctions
    [pt{1}, cd{1}] = ...
        Refine_TJunctions(pt{1}, TJ, cd2, ...
        curve, curve_num, curve_start, curve_end, curve_mode, EP);
    % pt : n by 2 matrix containing the positions of the detected corners,
    %      where n is the number of detected corners
    % cd : cpda curvature values of the detected corners
    
    if nargout >= 2
        varargout{1} = false(size(BW));
        varargout{1}(sub2ind([sizex,sizey],pt{1}(:,1),pt{1}(:,2))) = true;
        if nargout == 3,  varargout{2}{1} = cd;  end
    end
    
    
else
    pt{1} = [];
    for i=1:nargout-1, varargout{i} = [];   end
end

end % end of fastcpda_base


%% Subfunctions

%%
% |SELECTPOINT| - Select points from extracted curves.
%--------------------------------------------------------------------------
function [point_selected smoothed_curve Width Sig] = ...
    selectPoint(curve, curve_mode, curve_num, sizex, sizey)

s = 3.0; % minor axis
% S = 1.5*s; % in that case, 1.5 denotes the minimum ratio of major axis to
% minor axis of an ellipse, whose vertex could be detected as a corner
S = 4.0; % major axis
[gau w] = find_Gaussian(s);
[Gau W] = find_Gaussian(S);

extra = W-w;
gau1 = [zeros(1,extra) gau zeros(1,extra)];
DoG = Gau-gau1;
%t = 0.1;

smoothed_curve = cell(curve_num);
point_selected = cell(curve_num);

Width = w; Sig = s;
for i = 1:curve_num    
    x = curve{i}(:,2) - sizey/2;
    y = sizex/2 - curve{i}(:,1);
    L = size(x,1);
    if (L>W)
        % Calculate curvature
        if strcmpi(curve_mode(i,:),'loop')
            x1=[x(L-W+1:L);x;x(1:W)];
            y1=[y(L-W+1:L);y;y(1:W)];
        else
            x1=[ones(W,1)*2*x(1)-x(W+1:-1:2); ...
                x; ...
                ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
            y1=[ones(W,1)*2*y(1)-y(W+1:-1:2); ...
                y; ...
                ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
        end
        
        xx = conv(x1,DoG);
        xx = xx(W+1:L+3*W);
        yy = conv(y1,DoG);
        yy = yy(W+1:L+3*W);
        K = xx.^2 + yy.^2;
        
        % Find curvature local maxima as corner candidates
        N = size(K,1);
        n = 0;
        Search = 1;
        extremum = zeros(1,N);
        
        for j=1:N-1
            if (K(j+1)-K(j))*Search>0
                n=n+1;
                % in extremum, odd points is minima and even points is maxima
                extremum(n) = j;
                Search = -Search;
            end
        end
        if mod(n,2)==0
            n = n+1;
            extremum(n) = N;
        end
        extremum(n+1:N) = [];
        %n=size(extremum,2);
        %flag=ones(size(extremum));
        
        % Compare with adaptive local threshold to remove round corners
        %for j=2:2:n
        %    if K(extremum(j))<t
        %        flag(j)=0;
        %    end
        %end
        extremum = extremum(2:2:n);
        %flag=flag(2:2:n);
        %extremum=extremum(find(flag==1));
        extremum = extremum-W;
        extremum = extremum(extremum>0 & extremum<=L);
        
        xx = conv(x1,gau);
        xx = xx(W+1:L+3*W);
        yy = conv(y1,gau);
        yy = yy(W+1:L+3*W);
        
        smoothed_curve{i} = [xx yy];
        point_selected{i} = extremum;
        %Width = [Width W];
        %Sig = [Sig s];
    else
        smoothed_curve{i} = [];
        point_selected{i} = [];
        %Width = [Width 0];
        %Sig = [Sig 0];
    end
    %here = 1;
end
end % end of selectPoint


%%
% |GETCORNER| - Extract corners from curve representation.
%--------------------------------------------------------------------------
function [corners index cd] = ...
    getcorner(curve, curve_mode, curve_start, curve_num, T_angle, point_selected, sm_curve, W)

corners = [];
cd = [];

CLen = [10 20 30];
T = 0.2; % define the curvature threshold

index = cell(curve_num);

for i=1:curve_num;
    % C = [];  C3 = [];
    %x = curve{i}(:,2) - sizey/2;
    %y = sizex/2 - curve{i}(:,1);
    %curveLen = size(x,1);
    %[sig] = find_sig(curveLen);
    % smooth the curve with Gaussian kernel
    %[xs ys gau W] = smoothing(x,y,curveLen,curve_mode(i,:),sig,1);
    %xs = sm_curve{i}(:,1);
    %ys = sm_curve{i}(:,2);
    %W = Width(i);
    
    
    if ~isempty(sm_curve{i})
        xs = sm_curve{i}(:,2) ;
        ys = sm_curve{i}(:,1);
        L = size(xs,1);
        curveLen = L-2*W;
        %curveLen = L-2*W;
        %sig = Sigma(i);
        if L>1
            %if curve_mode(i,:)=='loop'
            %    xs1=[xs(curveLen-W+1:curveLen);xs;xs(1:W)];
            %    ys1=[ys(curveLen-W+1:curveLen);ys;ys(1:W)];
            %else %expand the ends to gaussian window
            %    xs1=[ones(W,1)*2*xs(1)-xs(W+1:-1:2); ...
            %         xs; ...
            %         ones(W,1)*2*xs(curveLen)-xs(curveLen-1:-1:curveLen-W)];
            %    ys1=[ones(W,1)*2*ys(1)-ys(W+1:-1:2); ...
            %         ys; ...
            %         ones(W,1)*2*ys(curveLen)-ys(curveLen-1:-1:curveLen-W)];
            %end
            %xs = xs1;
            %ys = ys1;
            %L = curveLen+2*W;
            extremum = point_selected{i};
            if size(extremum,2)>0
                C3 = zeros(3,L);
                for j = 1:3
                    chordLen = CLen(1,j);
                    C3(j,1:L) = abs(accumulate_chord_distance(xs,ys,chordLen,L,extremum,W));
                end
                c1 = C3(1,W+1:curveLen+W)/max(C3(1,W+1:curveLen+W));
                c2 = C3(2,W+1:curveLen+W)/max(C3(2,W+1:curveLen+W));
                c3 = C3(3,W+1:curveLen+W)/max(C3(3,W+1:curveLen+W));
                
                C = c1.*c2.*c3;
                %A = mean(C);
                L = curveLen;
                xs = xs(W+1:L+W);
                ys = ys(W+1:L+W);
                %flag = (extremum > W & extremum <= L+W);
                %extremum = extremum(flag == 1);
                %extremum = extremum-W;
                
                % Find curvature local maxima as corner candidates
                %extremum=[];
                %N=size(C,2);
                %n=0;
                %Search=1;
                
                %for j=1:N-1
                %    if (C(j+1)-C(j))*Search>0
                %        n=n+1;
                % % In extremum, odd points are minima and even points are maxima
                %        extremum(n)=j;
                % % minima: when K starts to go up; maxima: when K starts to go down
                %        Search=-Search;
                %    end
                %end
                %if mod(size(extremum,2),2)==0 %to make odd number of extrema
                %    n=n+1;
                %    extremum(n)=N;
                %end
                
                % accumulate candidate corners
                %n = size(extremum,2);
                %for j = 1:n
                %    cor = [cor; curve{i}(extremum(j),:)];
                %end
                
                n = size(extremum,2);
                flag = ones(size(extremum));
                
                % Compare each maxima with its contour average
                % if the maxima is less than local minima, remove it as false corner
                for j=1:n
                    if (C(extremum(j)) > T),   flag(j)=0;  end
                end
                %extremum = extremum(2:2:n); % only maxima are corners, not minima
                %flag = flag(2:2:n);
                extremum = extremum(flag==0);
                
                % Check corner angle to remove false corners due to boundary noise and trivial details
                %fl = 0;
                %if fl
                %flag=0;
                smoothed_curve = [xs,ys];
                while sum(flag==0) > 0
                    n = size(extremum,2);
                    flag = ones(size(extremum));
                    for j=1:n
                        % second argument of curve_tangent function is always
                        % the
                        % position of the extrema in the first argument which is
                        % an array of points between two extrema
                        if j==1 && j==n
                            ang = curve_tangent(smoothed_curve(1:L,:), extremum(j));
                        elseif j==1
                            ang = curve_tangent(smoothed_curve(1:extremum(j+1),:), extremum(j));
                        elseif j==n
                            ang = curve_tangent(smoothed_curve(extremum(j-1):L,:), ...
                                extremum(j)-extremum(j-1)+1);
                        else
                            ang = curve_tangent(smoothed_curve(extremum(j-1):extremum(j+1),:), ...
                                extremum(j)-extremum(j-1)+1);
                        end
                        % if angle is between T_angle = 162 and (360-T_angle) = 198
                        if ang>T_angle && ang<(360-T_angle)
                            flag(j) = 0;
                        end
                    end
                    
                    if size(extremum,2) == 0
                        extremum = [];
                    else
                        extremum = extremum(flag~=0);
                    end
                end
                % find corners which are not endpoints of the curve
                %extremum=extremum(find(extremum>0 & extremum<=curveLen));
                index{i} = extremum';
                % Sig(i,1) = sig;
                n = size(extremum,2);
                
                for j = 1:n
                    corners = [corners; curve{i}(extremum(j),:)];      %#ok
                    cd = [cd; C(extremum(j))];                         %#ok
                end
                                
                fl = 1;
                if fl && strcmp(curve_mode(i,:),'loop') && n>1
                    comp_corner = corners-ones(size(corners,1),1)*curve_start(i,:);
                    comp_corner = comp_corner.^2;
                    comp_corner = comp_corner(:,1) + comp_corner(:,2);
                    if min(comp_corner)>100
                        % add end points far from detected corners, i.e.
                        % outside of 5 by 5 neighbor
                        left = smoothed_curve(extremum(1):-1:1,:);
                        right = smoothed_curve(end:-1:extremum(end),:);
                        % detect corner at the first point or last point of the
                        % loop curve
                        ang = curve_tangent([left;right],extremum(1));
                        
                        % if angle is between T_angle = 162 and (360-T_angle) = 198
                        if ang>T_angle && ang<(360-T_angle)
                        else%if C(W+1)>T/2
                            corners = [corners; curve_start(i,:)];     %#ok
                            cd = [cd;5];                               %#ok
                        end
                    end
                end
            end
        end
    end
end

end % end of getcorner


%%
% |ACCUMULATE_CHORD_DISTANCE| - Accumulate chord distances.
%--------------------------------------------------------------------------
function Cd = ...
    accumulate_chord_distance(xs, ys, chordLen, curveLen, point_selected, W)

Cd = zeros(1,curveLen);

for j = 1:size(point_selected,2)
    k = point_selected(j)+W;
    %if k>W
    xk = xs(k); % (x1,y1) = point at which distance will be accumulated
    yk = ys(k);
    
    if k-chordLen+1 < 1,    s = 1;
    else                    s = k-chordLen+1;
    end
    
    for i = s:k-1
        if i+chordLen <= curveLen
            % (leftx,lefty) = current left point for which distance will be
            % accumulated
            x1 = xs(i);
            y1 = ys(i);
            
            % (rightx,righty) = current right point for which distance will
            % be accumulated
            x2 = xs(i+chordLen);
            y2 = ys(i+chordLen);
            
            % coefficients of st. line through points (x1,y1) and (x2,y2)
            a = y2-y1;
            b = x1-x2;
            c = x2*y1 - x1*y2;
            
            d = sqrt(a*a+b*b);
            if d~=0,  dist = (a*xk + b*yk + c)/d;
            else      dist = 0;
            end
            Cd(1,k) = Cd(1,k)+ dist;
        else
            break;
        end
    end
    
end

end % end of accumulate_chord_distance


%%
% |EXTRACT_CURVE| - Extract curves from input binary edge map: if the endpoint
% of a contour is nearly connected to another endpoint, fill the gap and continue
% the extraction.
%--------------------------------------------------------------------------
function [curve, curve_start, curve_end, curve_mode, cur_num, TJ, BW_edge] = ...
    extract_curve(BW, Gap_size)
% the default gap size is 1 pixel

[L,W] = size(BW);
BW1 = zeros(L+2*Gap_size,W+2*Gap_size);
BW_edge = zeros(L,W);
BW1(Gap_size+1:Gap_size+L,Gap_size+1:Gap_size+W) = BW;
[r,c] = find(BW1==1); % returns indices of non-zero elements
cur_num = 0;

while size(r,1)>0 % when number of rows > 0
    point = [r(1),c(1)];
    cur = point;
    % mask the pixel
    BW1(point(1),point(2)) = 0; 
    % find if any pixel around the current point is an edge pixel
    [I,J] = find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1); 
    while size(I,1) > 0 %if number of row > 0
        dist = (I-Gap_size-1).^2 + (J-Gap_size-1).^2;
        [~,index] = min(dist);
        p = point+[I(index),J(index)];
        % next is the current point
        point = p-Gap_size-1;
        % add point to curve 
        cur = [cur;point];                                             %#ok
        % mask the pixel
        BW1(point(1),point(2)) = 0;
        % find if any pixel around the current point is an edge pixel
        [I,J] = find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);              
    end
    
    % Extract edge towards another direction
    point  =[r(1),c(1)];
    BW1(point(1),point(2)) = 0;
    [I,J] = find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    
    while size(I,1)>0
        dist = (I-Gap_size-1).^2 + (J-Gap_size-1).^2;
        [~,index] = min(dist);
        point = point+[I(index),J(index)]-Gap_size-1;
        cur = [point;cur];                                             %#ok
        BW1(point(1),point(2)) = 0;
        [I,J] = find(BW1(point(1)-Gap_size:point(1)+Gap_size,point(2)-Gap_size:point(2)+Gap_size)==1);
    end
        
    % ?!!! for 512 by 512 image, choose curve if its length > 40 ?!!!
    if size(cur,1)>(size(BW,1)+size(BW,2))/25 
        % one can change this value to control the length of the extracted edges
        cur_num = cur_num+1;
        curve{cur_num} = cur-Gap_size;                                 %#ok
    end
    [r,c] = find(BW1==1);
    
end

curve_mode = char(zeros(cur_num,4));
for i=1:cur_num
    curve_start(i,:) = curve{i}(1,:);                                  %#ok
    curve_end(i,:) = curve{i}(size(curve{i},1),:);                     %#ok
    if (curve_start(i,1)-curve_end(i,1))^2+...
        (curve_start(i,2)-curve_end(i,2))^2<=25  %if curve's ends are within sqrt(32) pixels
        curve_mode(i,:) = 'loop';
    else
        curve_mode(i,:) = 'line';
    end
    BW_edge(curve{i}(:,1)+(curve{i}(:,2)-1)*L) = 1;
end

% if a contour goes just outsize of ends, i.e., outside of gapsize we note
% a T-junction there
if cur_num>0
    TJ = find_TJunctions(curve, cur_num, Gap_size+1); 
else
    curve{1} = [];
    curve_start = [];
    curve_end = [];
    curve_mode = [];
    cur_num = [];
    TJ = [];    
end

end % end of extract_curve


%%
% |FIND_TJUNCTIONS| - Find T-junctions in planar curves within (gap by gap)
% neighborhood, where gap = Gap_size + 1; edges were continued when ends are
% within (Gap_size by Gap_size)
%--------------------------------------------------------------------------
function TJ = find_TJunctions(curve, cur_num, gap)

TJ = [];

for i = 1:cur_num
    cur = curve{i};
    szi = size(cur,1);
    for j = 1:cur_num
        if i ~= j
            temp_cur = curve{j};
            comp_send = temp_cur - ones(size(temp_cur, 1),1)* cur(1,:);
            comp_send = comp_send.^2;
            comp_send = comp_send(:,1)+comp_send(:,2);
            % add curve strat-points as T-junctions using a (gap by gap)
            % neighborhood
            if min(comp_send)<=gap*gap
                TJ = [TJ; cur(1,:)];                                   %#ok
            end
            
            comp_eend = temp_cur - ones(size(temp_cur, 1),1)* cur(szi,:);
            comp_eend = comp_eend.^2;
            comp_eend = comp_eend(:,1)+comp_eend(:,2);
            % add end-points T-junctions using a (gap by gap) neighborhood
            if min(comp_eend) <= gap*gap
                TJ = [TJ; cur(szi,:)];                                 %#ok
            end
        end
    end
end
end % end of find_TJunctions


%%
% |REFINE_TJUNCTIONS| - Compare T-junctions with obtained corners and add 
% T-junctions to corners which are far away (outside a 5 by 5 neighborhood) 
% from detected corners
%--------------------------------------------------------------------------
function [corner_final c3] = ...
    Refine_TJunctions(corner_out, TJ, c2,curve, curve_num, curve_start, curve_end, curve_mode,EP)

% corner_final = corner_out;
c3 = c2;

% add end points
if EP
    corner_num = size(corner_out,1);
    for i=1:curve_num
        if size(curve{i},1)>0 && strcmpi(curve_mode(i,:),'line')
            
            % Start point compare with detected corners
            comp_corner = corner_out - ones(size(corner_out,1),1)*curve_start(i,:);
            comp_corner = comp_corner.^2;
            comp_corner = comp_corner(:,1) + comp_corner(:,2);
            if min(comp_corner) > 100       % Add end points far from detected corners
                corner_num = corner_num + 1;
                corner_out(corner_num,:) = curve_start(i,:);
                c3 = [c3;8];                                           %#ok
            end
            
            % End point compare with detected corners
            comp_corner = corner_out - ones(size(corner_out,1),1)*curve_end(i,:);
            comp_corner = comp_corner.^2;
            comp_corner = comp_corner(:,1) + comp_corner(:,2);
            if min(comp_corner) > 100
                corner_num = corner_num + 1;
                corner_out(corner_num,:) = curve_end(i,:);
                c3 = [c3;9];                                           %#ok
            end
        end
    end
end

% add T-junctions
corner_final = corner_out;

for i=1:size(TJ,1)
    % T-junctions compared with detected corners
    if size(corner_final)>0
        comp_corner = corner_final - ones(size(corner_final,1),1)*TJ(i,:);
        comp_corner = comp_corner.^2;
        comp_corner = comp_corner(:,1) + comp_corner(:,2);
        if min(comp_corner) > 100       % Add end points far from detected corners, i.e. outside of 5 by 5 neighbor
            corner_final = [corner_final; TJ(i,:)];                    %#ok
            c3 = [c3;10];                                              %#ok
        end
    else
        corner_final = [corner_final; TJ(i,:)];                        %#ok
        c3 = [c3;10];                                                  %#ok
    end
end
end % end of Refine_TJunctions
 

%%
% |CURVE_TANGENT| - Compute the tangent direction over a curve.
%--------------------------------------------------------------------------
function ang = curve_tangent(cur, center) 
% center is always the position of the corresponding extrema in cur

dir = zeros(2,1);
for i=1:2
    if i == 1,     curve = cur(center:-1:1,:);
    else           curve = cur(center:size(cur,1),:);  
    end
    L = size(curve,1);
    
    if L>3
        if sum(curve(1,:) ~= curve(L,:))~=0 % if not collinear
            M = ceil(L/2);
            x1 = curve(1,1);    y1 = curve(1,2);
            x2 = curve(M,1);    y2 = curve(M,2);
            x3 = curve(L,1);    y3 = curve(L,2);
        else
            M1 = ceil(L/3);     M2 = ceil(2*L/3);
            x1 = curve(1,1);    y1 = curve(1,2);
            x2 = curve(M1,1);   y2 = curve(M1,2);
            x3 = curve(M2,1);   y3 = curve(M2,2);
        end
        
        if abs((x1-x2)*(y1-y3) - (x1-x3)*(y1-y2))<1e-8  % straight line
            tangent_dir = ...
                angle(complex(curve(L,1)-curve(1,1), curve(L,2)-curve(1,2)));
        else
            % Fit a circle 
            x0 = 1/2 * (-y1*x2^2 + y3*x2^2 - y3*y1^2 - y3*x1^2 - y2*y3^2 + ...
                x3^2*y1 + y2*y1^2 - y2*x3^2 - y2^2*y1 + y2*x1^2 + y3^2*y1 + ...
                y2^2*y3) / (-y1*x2 + y1*x3 + y3*x2 + x1*y2 - x1*y3 - x3*y2);
            y0 = -1/2 * (x1^2*x2 - x1^2*x3 + y1^2*x2 - y1^2*x3 + x1*x3^2 - ...
                x1*x2^2 - x3^2*x2 - y3^2*x2 + x3*y2^2 + x1*y3^2 - x1*y2^2 + ...
                x3*x2^2) / (-y1*x2 + y1*x3 + y3*x2 + x1*y2 - x1*y3 - x3*y2);
            % R = (x0-x1)^2+(y0-y1)^2;

            radius_dir = angle(complex(x0-x1,y0-y1));
            if radius_dir<0
                radius_dir = 2*pi-abs(radius_dir);
            end
            
            adjacent_dir = angle(complex(x2-x1,y2-y1));
            
            if adjacent_dir<0
                adjacent_dir = 2*pi-abs(adjacent_dir);
            end
            
            tangent_dir = sign(sin(adjacent_dir-radius_dir))*pi/2 + radius_dir;
            if tangent_dir<0
                tangent_dir = 2*pi-abs(tangent_dir);
            elseif tangent_dir>2*pi
                tangent_dir = tangent_dir-2*pi;
            end
        end
    
    else % very short line
        tangent_dir = ...
            angle(complex(curve(L,1)-curve(1,1),curve(L,2)-curve(1,2)));
    end
    dir(i) = tangent_dir*180/pi;
end

ang = abs(dir(1) - dir(2));

end % end of curve_tangent


%%
% |MAKEGFILTER|
%--------------------------------------------------------------------------
function [G W] = makeGFilter(sig)

GaussianDieOff = .0001; 
pw = 1:100;

ssq = sig*sig;
W = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff, 1, 'last');
if isempty(W),     W = 1;   end

t = (-W:W);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq); 
G = gau/sum(gau);

end % end of makeGFilter


%%
% |FIND_GAUSSIAN|
%--------------------------------------------------------------------------
function [gau width] = find_Gaussian(sig)
GaussianDieOff = .0001; 
pw = 1:30; 
ssq = sig*sig;

width = find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff, 1, 'last');
if isempty(width),    width = 1;   end

t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq); 
gau = gau/sum(gau);

end % end of find_Gaussian


%%
% |ENLARGE|
%--------------------------------------------------------------------------
function [xse yse] = enlarge(xs, ys, CL, curve_mode)                   %#ok
% CL = chord length
L = size(xs,1);

if strcmpi(curve_mode,'loop') % wrap around the curve by CL pixles at both ends
    xse = [xs(L-CL+1:L); ...
        xs; ...
        xs(1:CL)];
    yse = [ys(L-CL+1:L); ...
        ys; ...
        ys(1:CL)];
else % extend each line curve by CL pixels at both ends
    xse = [ones(CL,1)*2*xs(1)-xs(CL+1:-1:2); ...
        xs; ...
        ones(CL,1)*2*xs(L)-xs(L-1:-1:L-CL)];
    yse = [ones(CL,1)*2*ys(1)-ys(CL+1:-1:2); ...
        ys; ...
        ones(CL,1)*2*ys(L)-ys(L-1:-1:L-CL)];
end
end % end of enlarge


%%
% |SMOOTHING|
%--------------------------------------------------------------------------
function [xs ys gau W] = smoothing(x, y, L, curve_mode, sig, mode)     %#ok
[gau W] = makeGFilter(sig);

if L>W
    if strcmpi(curve_mode,'loop') % wrap around the curve by W pixles at both ends
        x1 = [x(L-W+1:L); ...
            x; ...
            x(1:W)];
        y1 = [y(L-W+1:L); ...
            y; ...
            y(1:W)];
    else % extend each line curve by W pixels at both ends
        x1 = [ones(W,1)*2*x(1)-x(W+1:-1:2); ...
            x; ...
            ones(W,1)*2*x(L)-x(L-1:-1:L-W)];
        y1 = [ones(W,1)*2*y(1)-y(W+1:-1:2); ...
            y; ...
            ones(W,1)*2*y(L)-y(L-1:-1:L-W)];
    end
    
    xx = conv(x1,gau);    
    yy = conv(y1,gau);
    if (mode == 1)
        xs = xx(W+1:L+3*W);
        ys = yy(W+1:L+3*W);    
    else
        xs = xx(2*W+1:L+2*W);
        ys = yy(2*W+1:L+2*W);    
    end
else
    xs = [];
    ys = [];    
end
end % end of smoothing


%%
% |FIND_SIG|
%--------------------------------------------------------------------------
function [sig] = find_sig(L)                                           %#ok
if L<=100
    sig = 3;
    %wid = 4;
elseif L<=200
    sig = 3;
    %wid = 8;
else
    sig = 3;
    %wid = 12;
end
end % end of find_sig


