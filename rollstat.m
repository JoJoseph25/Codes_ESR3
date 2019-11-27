%	Differs from "most" rolling window tools in that this function operates
%	on all x-values within the specified window, rather than just rolling
%	by a specified number of x-index values.
% Copyleft :-)   2015 by Carl Witthoft  (carl@witthoft.com)

function [ rollstats] = rollstat( y, window, slide,funin)
% rollstat
% to run rolling window calculations.
% Inputs:
%	y - vector of your data, or nX2 (or 2Xn) array of y,x pairs.  
%       If x-values are included, they must be monotonic increasing.
%		User-supplied x-values can be a non-uniform time or distance sequence.
%		The window and slide sizes are adjusted to the 'true' time
%		or distance window size desired. The default x-vector 
%		contains equally spaced values.
%	window - how wide (how many samples or x-width) to crunch each time
%	slide - how many samples or x-width to move the window each time 
%	funin - optional: either a char string name of function or a function handle
%		to execute in place of the builtin default function. The function
%		must return a single item, be it scalar, struct, or cell (not sure
%		about cell).
% Output:
%	Internal default: A structure containing 
%   max,min,mean,median,sigma,skew for all	windows
%	When the user supplies a function, the output is a struct containing 
%   whatever that function generates, calculated at each window position.
%   In either case, the following info is returned as well:
%	-- 'actual' x-values for each window
%	-- the actual window and "slide" values
%	-- the indices at which the windows begin
%
%	All users are more than welcome to replace the provided default
%	function with a function of their choosing.
 
if ~exist('window','var') || (numel(window)< 1) , 
	window = floor(numel(y)/10); 
end;
if ~exist('slide','var') || (numel(slide)<1) , slide = 1; end;
switch (size(y,1)),
	case 1,
		x = 1:numel(y);
	case 2
		x = y(2,:);
		y = y(1,:);
	otherwise,
		switch (size(y,2)),
			case 1
				x = 1:numel(y);
			case 2
				x = y(:,2);
				y = y(:,1);
			otherwise
				error('Max of either 2 rows (y,x) or two columns (y,x) allowed');
		end
end
%next, process the input function name if there is one.
if exist('funin','var'),
	if ischar(funin);
		thefunc = eval(['@(x) ' funin '(x)']);		
	else
		if isa(funin,'function_handle'),
			thefunc = funin;
		else 
		error('"funin" must be function name or handle.');
		end
	end	
else
	thefunc = @(x)defaultfunc(x);
end

foowin = [];
fooslide = []; 
jj = 1;
kk = 1;
refidy(1) = 1; 
% limit processing window to actual extent of data
maxjj = find( (x(end)-x) >= window,1,'last'); 
% if there are dupes, this line picks them all up:
maxjj = find(x <= x(maxjj),1,'last');
 
while jj <= maxjj , 
	tfoo = find(x >= (x(jj)+window) ,1,'first');
	if numel(tfoo),
		foowin(kk)= tfoo - jj;
	else
		break
	end;
	tbar =  find(x >= (x(jj)+ slide),1,'first');
	if numel(tbar),  
		fooslide(kk) = tbar -jj ;
		jj = jj + fooslide(kk);
		kk = kk+1;
		refidy(kk) = refidy(kk-1)+fooslide(kk-1);
	else 
		break  
	end
end  %  end of  while loop
% remove last refidy since it's not real 
 refidy = refidy(1:(end-1));
% Notice that the next line defines the class of tstats
tstats = thefunc(NaN);
tstats = repmat(tstats, numel(refidy),1); 

% starting index of y is the sum of all trueslides from 1 to j-1
%  the first window starts at index 1
for j = 1 : numel(refidy),
	tstats(j) = thefunc( y(refidy(j) : ...
		( min(refidy(j)+foowin(j),numel(y)) ) )  ) ;	
end

rollstats.xout = x(refidy);
rollstats.stats=tstats;
rollstats.twin = foowin; 
rollstats.tslide = fooslide; 
rollstats.idx = refidy;
end

function astat = defaultfunc(x)
% semi-port of R-function mystat (CRAN package cgwtools)
x = reshape(x,1,numel(x));
x = x(~isnan(x));
skewx = theskew(x);
% A fix for when I enter just NaN or some other foolish thing.
if(~numel(x)), 
	xM=NaN;
	xm = NaN;
else
	xM = max(x);
	xm=min(x);
	end;
astat = struct('max',xM, 'min', xm, ...
	'mean',mean(x), 'median',median(x), 'sd', std(x), 'skew', skewx) ;
end
	
function skw = theskew(x) 
%# skew  from  the R CRAN e1071 package.
% use as internal func only
 skw = sum((x-mean(x)).^3)/(numel(x).*std(x).^3);
 end