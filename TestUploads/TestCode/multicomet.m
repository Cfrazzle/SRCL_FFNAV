function multicomet(varargin)
%MULICOMET  Extends COMET to display multiple trajectories. No tail.
%   MULTICOMET(Y) displays animated comet plots of the columns vectors of Y.
%   MULTICOMET(X,Y) displays animated comet plots of the column vectors of Y vs. X.
%
%   MULTICOMET(AX,...) plots into AX instead of GCA.
%
%   Example:
%   
%   multicomet([rand(100,1) + 5,rand(100,1), rand(100,1) - 5])
%
%
%   See also COMET, COMET3.


% Parse possible Axes input
[ax,args,nargs] = axescheck(varargin{:});

error(nargchk(1,3,nargs,'struct'));

% Parse the rest of the inputs
if nargs < 2, x = args{1}; y = x; x = repmat((1:size(y,1))',1,size(y,2)); end
if nargs == 2
    [x,y] = deal(args{:});
    if size(x,2) == 1
        x = repmat(x,1,size(y,2));
    end
end
p = .1;

ax = newplot(ax);
if ~ishold(ax)
    [minx,maxx] = minmax(x);
    [miny,maxy] = minmax(y);
    axis(ax,[minx maxx miny maxy])
end
co = get(ax,'colororder');

[m,n] = size(x);

for i=1:n
    % Choose first three colors for head, body, and tail
    head(i) = line('parent',ax,'color',co(1,:),'marker','o','erase','xor', ...
        'xdata',x(1,i),'ydata',y(1,i));
    body(i) = line('parent',ax,'color',co(i,:),'linestyle','-','erase','none', ...
        'xdata',[],'ydata',[]);
end

k = round(p*m);

% This try/catch block allows the user to close the figure gracefully
% during the comet animation.
try
    % Grow the body
    for i = 2:k+1
        j = i-1:i;
        for a = 1:n
            set(head(a),'xdata',x(i,a),'ydata',y(i,a))
            set(body(a),'xdata',x(j,a),'ydata',y(j,a))
        end
        drawnow
    end
    
    % Primary loop
    for i = k+2:m
        j = i-1:i;
        for a = 1:n
            set(head(a),'xdata',x(i,a),'ydata',y(i,a))
            set(body(a),'xdata',x(j,a),'ydata',y(j,a))
        end
        drawnow
    end

catch E
    if ~strcmp(E.identifier, 'MATLAB:class:InvalidHandle')
        rethrow(E);
    end
end

function [minx,maxx] = minmax(x)
minx = min(x(isfinite(x)));
maxx = max(x(isfinite(x)));
if minx == maxx
    minx = maxx-1;
    maxx = maxx+1;
end
