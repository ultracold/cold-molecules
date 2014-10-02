%% TEMPERATURE MEASUREMENT OF COLD MOLECULES
%  AUTHOR: ultracold
%  MATLAB VERSION: R2014a
%%-------------------------------------------------------------------------
%
%%-------------------------------------------------------------------------
%% CONTENTS
%
% ESTIMATING THE TEMPERATURE OF COLD MOLECULES
%   - SIMPLE DATA IMPORT
%   - PLOTTING DATA AND FOR LOOPS IN MATLAB
%   - NaN AND DATA
%   - FITTING DATA AND USING FUNCTIONS IN MATLAB
%   - USING FIT RESULTS WISELY WITH RESPECT TO EXPERIMENTAL DATA
%   - ESTIMATING A TEMPERATURE



%%-------------------------------------------------------------------------
%% 1. CLEAR WORKSPACE
% Close all figures, clear workspace, clear command window.
close all
clear all
clc



%%-------------------------------------------------------------------------
%% 2. FIGURE COUNT
% We will be plotting a few figures, so I will here create a figure counter
% variable to keep track of them.
figureCount = 0;



%%-------------------------------------------------------------------------
%% 5. ESTIMATING THE TEMPERATURE OF COLD MOLECULES
%
% We will now move on to working with some data taken from experiments
% concerning cold molecules trapped in a microchip environment.
% The molecules are captured in an array of microtraps that are spaced by
% 120 microns from one another. The molecules are then released from the
% microtraps and their free ballistic expansion is monitored over time.
% Through this expansion, it is possible to estimate the temperature of the
% cold molecular ensemble, which is what we will do here.
% The data analysed here are line profiles taken from integrated images of
% the cold molecules for different ballistic expansion times.
% The data are in a .txt file and are pairs of
% columns showing the position and the molecule count at that position.


%% SIMPLE DATA IMPORT
%
% We will begin by importing the data:

filename = 'ColdMoleculeData.txt';
molData = importdata(filename);

% Note: one should be careful with the function 'importdata', since,
% firstly, it is less advaced in older versions of MATLAB and,
% secondly, you have to be sure of what it is doing with the data it
% imports!

% From experiments, the ballistic expansion times corresponding to the data
% we just imported are:

expansionTimes = [9; 15; 19; 22]; % in microseconds

% Note: In total, there are 4 expansionTimes and 8 columns of data in
% molData, two for each expansionTime.


%% PLOTTING DATA AND FOR LOOPS IN MATLAB

% Create a figure window to plot our data in:
figureCount = figureCount +1;
figHandle = figure(figureCount);
set(figHandle,'name', 'Free Expansion of Cold Molecules', ...
    'Position', [100 300 1200 300])


% Here is a FOR loop in MATLAB.  It creates 4 subfigures in the figure
% window we just defined above using the molecule data we imported.
% It plots the position on the x-axis and the number of molecules at that
% position on the y-axis. 
for iSubFig = 1:4
    subplot(1, 4, iSubFig)
    plot(molData(:, iSubFig*2-1)*10^3, molData(:, iSubFig*2))
    axis square
    ylim([0 Inf])
    xlabel('Position in microns')
    ylabel('No. of Molecules')
    title(['Expansion Time = ' ...
          num2str(expansionTimes(iSubFig)) ...
          ' microsec'])
    hold on  % We will plot again in a few moments on the same figure. 
end


%% NaN and DATA
% When MATLAB imports data, it fills the empty items with NaN,
% since some data columns are longer than others and all are imported
% together. This command just takes away any NaN values and makes them
% zero.  In general, when working with data, setting NaN to zero is not
% always the best option, but in this case it suffices.
molData(isnan(molData)) = 0;


%% FITTING DATA AND USING FUNCTIONS IN MATLAB

% Here we will fit the experimental data (i.e. the line profiles of the
% cold molecular clouds) with a multi-Gaussian fit, i.e. we will
% assume that the molecular clouds are guassian in shape (which in this
% case is a good assumption). A multi-Gaussian in this case just means
% the fitting of the sum of several gaussian profiles at the same time.

% Here are the starting values for the cloud widths (standard deviation) 
% for use in the fit.
peakWidthGuess = [20, 27, 35, 40];


% Choose position of first peak to fit (which in this case is to the
% left of the data shown in the figure, so that we definitely
% capture all the peaks in the experimental data).
firstPeakPosition = [76600, 76640, 76580, 76660];

% Pre-allocation of a data array to store the fits to the data.
% Here I have chosen a cell array because this allows the storage
% of more complicated objects in each unit of the cell (for example,
% the cfit object created by a fit in MATLAB).
fitMolData = cell(1,4);         

% Perform the fit.
% Since we will fit the data 4 times, it makes sense to create
% a function that we can reuse for each data set, and to put
% this function in a separate .m file for code clarity.
% Here we loop through the 4 data sets, fitting to each.
for iSubFig = 1:4
     fitMolData{iSubFig} = multiGaussFit(molData(:, iSubFig*2-1)*10^3, ...
                           molData(:, iSubFig*2), ...
                           firstPeakPosition(iSubFig), ...
                           peakWidthGuess(iSubFig));
end         

% Now we can plot each of our fits on top of the experimental data
% we already displayed.
for iSubFig = 1:4
    subplot(1, 4, iSubFig)
    plot(fitMolData{iSubFig})
    xlim([(firstPeakPosition(iSubFig) + 120/2) ...
           firstPeakPosition(iSubFig) + (5.5)*120])
    ylim([0 Inf])
    xlabel('Position in microns')
    ylabel('No. of molecules')
    hold off 
end 


%% USING FIT RESULTS WISELY WITH RESPECT TO EXPERIMENTAL DATA

% Collect cloud width values from all fits.

% Initialize cell array to store fit widths of the molecule clouds.
fitWidths = cell(1,length(expansionTimes));

% Put fitted cloud widths into fitWidths cell array.
% (The multiGaussFit function used above fits 8 gaussians in total.)
for iSubFig = 1:4
    fitWidths{iSubFig} = [fitMolData{iSubFig}.aaa, ...
                          fitMolData{iSubFig}.bbb, ...
                          fitMolData{iSubFig}.ccc, ...
                          fitMolData{iSubFig}.ddd, ...
                          fitMolData{iSubFig}.eee, ...
                          fitMolData{iSubFig}.fff, ...
                          fitMolData{iSubFig}.ggg, ...
                          fitMolData{iSubFig}.hhh];
end

% Just using functions blindly on experimental data is never a
% good policy, as is blindly taking the output values as the 
% correct values. From our experimental work, we know that there
% is only a specific area of the images that is unaffected by
% (ion-) optical aberrations. This area includes the peaks 2:6 in
% the images (you will just have to trust me on this!)

% Pick widths of the clouds in the aberration-free are of the image:
peakChoice = {[2:6], [2:6], [2:6], [2:6]};

% Select the aberration-free cloud widths from fitWidths.
fitWidthsChoice = cell(1,4);
for iSubFig = 1:4
     fitWidthsChoice{iSubFig} = fitWidths{iSubFig}(peakChoice{iSubFig});
end     


%% ESTIMATING A TEMPERATURE

% The ballistic expansion of the cloud is (ideally) governed by
% the following relationship:
% cloudWidth^2 = initial_cloudWidth^2 + kB*T/mass * expansionTime^2
% where T is the temperature (kB is the Boltzmann constant and mass is mass
% is the mass of a molecule. Hence, the gradient of the relationship
% between the width squared and the expansion time squared gives the
% temperature, which we will estimate below.

% Here we calculate the mean square cloud width for each expansion time.
% We also calculate the std. dev. of the mean squared width, and hence
% the standard error.
meanSqWidth = zeros(4,1);
stdSqWidth = zeros(4,1);
stdErrorWidth = zeros(4,1);
for iSubFig = 1:4
    meanSqWidth(iSubFig) = mean(fitWidthsChoice{iSubFig}.^2);
    stdSqWidth(iSubFig) = std(fitWidthsChoice{iSubFig}.^2);
    stdErrorWidth(iSubFig) = stdSqWidth(iSubFig)./ ...
                             sqrt(length(fitWidthsChoice{iSubFig}));
end

% Now we will plot the mean square width against expansion time squared.
figureCount = figureCount +1;
figure(figureCount)
figHandle = figure(figureCount);
set(figHandle,'name', 'Expansion - Width squared against time squared')
errorbar(expansionTimes.^2,meanSqWidth,stdErrorWidth,'o')
xlabel('expansion time squared')
ylabel('mean square width')
title('Checking Linearity of Ballistic Expansion')
hold on  % We will plot on this figure again in a moment. 


% We want to fit a linear line to our data, which we do using a new
% function: linearFit.
expansionFit = linearFit(expansionTimes.^2, meanSqWidth);
plot(expansionFit)  % plot fit on figure.
hold on

% As can be seen in the figure, the fourth data point (i.e. for an 
% expansion time of 22 microseconds) has a relatively large standard error
% and the other three data points appear to be more aligned. Indeed,
% looking at Figure 4: Free Expansion of Cold Molecules, the data for the
% fourth expansion time is much noisier than for the other three expansion
% times (this is because the molecule clouds have spread out and the
% signal-to-noise has reduced).
% Therefore, now we will fit only to the first three expansion times:
expansionFitFirst3Times = linearFit(expansionTimes(1:3).^2, ...
                                    meanSqWidth(1:3));
plot(expansionFitFirst3Times, 'k')  % plot fit on figure.
xlabel('expansion time squared')
ylabel('mean square width')
title('Checking Linearity of Ballistic Expansion')
hold off

% As can be seen in the figure (Figure 5), the new linear fit to only the
% first three data points demonstrates the (expected) linear behaviour of
% the expansion.  Of course, it is open to interpretation as to whether the
% fourth data point is simply noisy (and therefore difficult to fit) or 
% whether something happens experimentally that is as yet unexplained.

% Either way, we will now convert both fits into temperatures:

% Here are some physical constants we will need.
massCO = 29*1.67*10^-27;  % mass of a carbon monoxide molecule.
kB = 1.38*10^-23;         % Boltzmann constant.
milliKelvin = 10^3;       % conversion factor from K to mK.
unitsConversion = massCO/kB*10^3;  % to convert fit values to temperatures.

% We can now calculate a temperature along with the corresponding errors,
% which we take from the fit using 'confint'.
% Firstly, we will do this for when fitting to all four expansion times:
linearFitErrors = confint(expansionFit);
temperature = expansionFit.a*unitsConversion;
temperatureErrorDown = linearFitErrors(1)*unitsConversion - temperature;
temperatureErrorUp = linearFitErrors(2)*unitsConversion - temperature;

% and then print the results to the command window.
fprintf(['\nThe temperature of the cold molecular cloud, \n(when ' ...
         'fitting all four expansion times) is\n' ...
         num2str(round(temperature)) ...
         ' mK \nwith errors of \n+' ...
         num2str(round(temperatureErrorUp)) ...
         ' mK \nand\n' ...
         num2str(round(temperatureErrorDown)) ...
         ' mK.\n'])

% Secondly, we will calculate the temperature with errors when fitting only
% to the first three expansion times:
linearFitErrorsFirst3Times = confint(expansionFitFirst3Times);
temperatureFirst3Times = expansionFitFirst3Times.a*unitsConversion;
temperatureErrorDownFirst3Times = ...
                     linearFitErrorsFirst3Times(1)*unitsConversion - ...
                     temperatureFirst3Times;
temperatureErrorUpFirst3Times = ...
                     linearFitErrorsFirst3Times(2)*unitsConversion - ...
                     temperatureFirst3Times;
                 
% and then print the results to the command window.
fprintf(['\nThe temperature of the cold molecular cloud, \n(when ' ...
         'fitting only the first three expansion times) is\n' ...
         num2str(round(temperatureFirst3Times)) ...
         ' mK \nwith errors of \n+' ...
         num2str(round(temperatureErrorUpFirst3Times)) ...
         ' mK \nand\n' ...
         num2str(round(temperatureErrorDownFirst3Times)) ...
         ' mK.\n'])




%% ------------------------------------------------------------------------
% end of script.

