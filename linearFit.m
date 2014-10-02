function [fitResult] = linearFit(x, y)

% linearfit is to fit a linear line y = a*x + b to x and y inputs.

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[0, -Inf],...
               'Upper',[Inf, Inf],...
               'StartPoint',[1 0]);
ft = fittype('a*x + c','options',fo); % linear function.

[fitResult] = fit(x, y, ft); % return fit result.
