function r2 = fit_sine(x,y,varargin)
% taken from https://www.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function

if length(varargin) >= 2
    per = varargin{2};
else
    per = 180;
end

yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
% per = 2*mean(diff(zx));                     % Estimate period
ym = mean(y);                               % Estimate offset
% fit = @(b,x)  b(1).*(sin(2*pi*x./b(2) + 2*pi/b(3))) + b(4);
fit = @(b,x)  b(1).*(sin(2*pi*x./per + 2*pi/b(2))) + b(3); % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
s = fminsearch(fcn, [yr; -1;  ym]);                       % Minimise Least-Squares
% xp = linspace(min(x),max(x));

yp = fit(s,x);

if length(varargin)>=1
    if varargin{1}
%         figure(1)
        plot(x,y,'b',  x,yp, 'r')
        grid
    end
end

r2 = 1-(sum(sqrt((y-yp).^2))/sum(sqrt((y-ym).^2)));