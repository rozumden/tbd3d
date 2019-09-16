function [idx1 idx2 img] = renderLineMulti(pts, sz)
% (POZN: mozna chyba, ze idx1/2 na vystupu nerespektuje sz, takze to vraci vzdy absolutni souradnice ktere nemusi korespondovat s vyrenderovanym obrazkem)
% Renders line or sequence of lines as image. Simple, no antialiasing.
%
% pts - (N,2) sequence of points (at least two) specified as [idx1 idx2] (not (x,y)) of end points (breaking points) of the line.
% sz - optional, if specified, also returs rendered line as image of size sz (otherwise smallest size is used and 'pts' then has only relative meaning)
%
% idx1/2 - row/col indices of nonzero line pixels, in order.
% img - if specified, also returns line as logical image

idx1 = []; idx2 = [];
for i=1:size(pts,1)-1
    [i1 i2] = renderLineSingle(pts(i,:), pts(i+1,:));
    if(i > 1)
        idx1 = [idx1; i1(2:end)]; idx2 = [idx2; i2(2:end)]; % do not repeat line start
    else
        idx1 = [idx1; i1]; idx2 = [idx2; i2];
    end
end

% render image if required
if(nargout >= 3)
    if(nargin >= 2 && ~isempty(sz))
        offset = [0 0];
    else
        offset = [min(idx1) min(idx2)]-1;
        sz = [max(idx1) max(idx2)]-offset;
    end

    img = false(sz);
	idx = [idx1 idx2]-offset;
	idx = idx(all(idx > 1 & idx <= sz,2),:);
    %idx = sub2ind(sz, idx1-offset(1), idx2-offset(2));
    img((idx(:,2)-1)*sz(1)+idx(:,1)) = true;
end
end


function [idx1 idx2] = renderLineSingle(from, to)
% line drawing algorithm from file exchange (https://www.mathworks.com/matlabcentral/fileexchange/28190-bresenham-optimized-for-matlab)
% no antialiasing, returns indices of line pixels [idx1(i) idx2(i)] is i-th line pixel
x1=round(from(2)); x2=round(to(2));
y1=round(from(1)); y2=round(to(1));
dx=abs(x2-x1);
dy=abs(y2-y1);
steep=abs(dy)>abs(dx);
if steep t=dx;dx=dy;dy=t; end

%The main algorithm goes here.
if dy==0 
    q=zeros(dx+1,1);
else
    q=[0;diff(mod([floor(dx/2):-dy:-dy*dx+floor(dx/2)]',dx))>=0];
end
%and ends here.

if steep
    if y1<=y2 idx1=[y1:y2]'; else idx1=[y1:-1:y2]'; end
    if x1<=x2 idx2=x1+cumsum(q);else idx2=x1-cumsum(q); end
else
    if x1<=x2 idx2=[x1:x2]'; else idx2=[x1:-1:x2]'; end
    if y1<=y2 idx1=y1+cumsum(q);else idx1=y1-cumsum(q); end
end
end