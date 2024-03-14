function strataWorkbench
% Code for StrataWorkbench
% Written by P Burgess, University of Liverpool
% Started 2/4/2016, incoporating various older code into a new GUI and creating basis for developing new functionality

    clear all;
    
    gui.main = 0;
    gui.f1 = 0;

    data.loaded = 0;    % Boolean flag to show if data loaded or not; if not, subsequent functions disabled
    data.sectionName = '';
    data.faciesCodes = zeros(1,1);
    data.faciesThick = zeros(1,1);
    data.sectLength = 0.0;
    data.maxNumbOfFacies = 0;
    data.sectionMean = 0.0;
    
    results.markovOrderMetric = zeros(1,1);

    initializeGUI(gui, data, results);
end

function initializeGUI(gui, data, results)

 
    % ScreenSize is a four-element vector: [left, bottom, width, height]:
    scrsz = get(0,'ScreenSize'); % vector 
    
    
    %% Create the main graphics window, with the vertical section plot, a selection of other output, and the main buttons and text parameter inputs

    scrWidthProportion = 0.75;
    maxButtonCount = 15; % Maximum number of buttons and text boxes we want to fit down right hand side of the screen
    scrHeightIncrement = (scrsz(4)*0.75) / maxButtonCount; % Use this to space controls down right side of the main window
    controlStartY = (scrsz(4) * 0.8) - (scrHeightIncrement / 2);
    controlStartX = (scrsz(3) * scrWidthProportion) - 420;
    
    % position requires left bottom width height values. screensize vector
    % is in this format 1=left 2=bottom 3=width 4=height
    % Hide the window as it is constructed by setting visible to off
    gui.main = figure('Visible','off','Position',[1 scrsz(4)*scrWidthProportion scrsz(3)*scrWidthProportion scrsz(4)*0.8]);
    
       %  Construct the control panel components.
       buttonCount = 1; % Count the elements as they are created, to control the position they appear on the screen
       
       uicontrol('style','text','string','Data folder path:','Position',[controlStartX+40, controlStartY-(scrHeightIncrement*buttonCount), 200, 15]);
       hSectionFpath = uicontrol('Style','edit','String','sectionData/','Position',[controlStartX+200, controlStartY-(scrHeightIncrement*buttonCount), 200, 25]);
       buttonCount = buttonCount + 1;
       
       uicontrol('style','text','string','Section data filename:','Position',[controlStartX+40, (controlStartY-scrHeightIncrement*buttonCount), 200, 15]);
       hSectionFname = uicontrol('Style','edit','String','exampleLog1.txt','Position',[controlStartX+200, (controlStartY-scrHeightIncrement*buttonCount), 200, 25]);
       buttonCount = buttonCount + 1;
       
       uicontrol('style','text','string','Section colour codes filename:','Position',[controlStartX+20, controlStartY-(scrHeightIncrement*buttonCount), 200, 15]);
       hSectionColourMapFname = uicontrol('Style','edit','String','exampleLog1lithoInfo.txt','Position',[controlStartX+200, controlStartY-(scrHeightIncrement*buttonCount), 200, 25]);
       buttonCount = buttonCount + 1;
       
       uicontrol('style','text','string','Bed count marker:','Position',[controlStartX+40, controlStartY-(scrHeightIncrement*buttonCount), 200, 15]);
       hSectionBedCountData = uicontrol('Style','edit','String','10','Position',[controlStartX+200, controlStartY-(scrHeightIncrement*buttonCount), 200, 25]);
       buttonCount = buttonCount + 1;
   
   uicontrol('Style','pushbutton','String','Load and plot section data',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.6 1.0 0.6],...
          'Callback',{@loadSectionButton_callback});
      buttonCount = buttonCount + 1;
      
   uicontrol('Style','pushbutton','String','Plot section data in scroll window',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.6 1.0 0.6],...
          'Callback',{@plotSectionWithScrollButton_callback})
      buttonCount = buttonCount + 1;
      
   uicontrol('Style','pushbutton','String','Bed thickness bar chart',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.6 1.0 0.6],...
          'Callback',{@plotBedThicknessBarChart_callback})
      buttonCount = buttonCount + 1;
      
    uicontrol('Style','pushbutton','String','Plot TP Matrix for observed strata',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@plotTPMatrix_callback});
      buttonCount = buttonCount + 1;
      
    uicontrol('Style','pushbutton','String','Plot facies frequency histogram',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@plotFaciesFrequencyHistogram_callback});
      buttonCount = buttonCount + 1;
      
    uicontrol('Style','pushbutton','String','Calculate and plot random models',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@CalculateAndPlotRandomModel_callback});
      buttonCount = buttonCount + 1;
      
    uicontrol('Style','pushbutton','String','Calculate and plot ideal cycles',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@CalculateAndPlotOptimisedCycles_callback});
      buttonCount = buttonCount + 1;
      
    uicontrol('Style','pushbutton','String','Moving window analysis',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@movingWindowAnalysisShowWindow_callback});
      buttonCount = buttonCount + 1;
      
    uicontrol('Style','pushbutton','String','Spectral analysis',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*buttonCount),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@spectralAnalysisShowWindow_callback});
     
      
    uicontrol('Style','pushbutton','String','Close any windows and clear data',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*maxButtonCount),200,25],...
          'BackgroundColor',[1.0 0.4 0.4],...
          'Callback',{@resetButton_callback});
    
   % Assign the GUI a name to appear in the window title.
   set(gui.main,'Name','Strata Workbench')
   % Move the GUI to the center of the screen.
   movegui(gui.main,'center')
   % Make the GUI visible.
   set(gui.main,'Visible','on');
   
   gui.movingWindowAnalysisWindowVisible = 0;

    function loadSectionButton_callback(source, eventdata)
        fprintf('\n\n===================================================== Loading Data =====================================================\n');
        
        dataFileNameAndPath = strcat(get(hSectionFpath,'String'),get(hSectionFname,'String'));
        coloursFileNameAndPath = strcat(get(hSectionFpath,'String'),get(hSectionColourMapFname,'String'));
        bedCountMarker = str2double(get(hSectionBedCountData,'string')); 

        [data, success] = loadSectionData(data, dataFileNameAndPath, coloursFileNameAndPath);
        if success == 1
            gui = plotObservedSectionInMainWindow(gui, data, bedCountMarker);
            data.loaded = 1; % Set flag to enable rest of the button functions
        end
    end

    function plotSectionWithScrollButton_callback(source, eventdata)
        fprintf('\n\n===================================================== Plotting Scrollable Section Data =====================================================\n');
        
        bedCountMarker = str2double(get(hSectionBedCountData,'string')); 
        
        if data.loaded
            gui = plotObservedSectionInScrollWindow(gui, data, bedCountMarker);
        end
    end

    function plotBedThicknessBarChart_callback(source, eventdata)
        fprintf('\n\n===================================================== Plotting Section Data as Bed Thickness bar chart =====================================================\n');
        
        bedCountMarker = str2double(get(hSectionBedCountData,'string')); 
        
        if data.loaded
            gui = plotBedThicknessBarChart(gui, data, bedCountMarker);
        end
    end

    function plotTPMatrix_callback(source,eventdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Plot TP Matrix for Observed Strata =====================================================\n');
            gui.sp2 = subplot('Position',[0.2 0.55 0.3 0.4]);
            [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 1, gui.sp2); % 1 is plot flag, so plot the results

            fprintf('Section gives m statistic %5.4f and one-offset diagonal check %5.4f\n', markovOrderMetric, oneOffsetMValueDiagCheck);
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function plotFaciesFrequencyHistogram_callback(source,eventdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Plot Facies Frequency Histogram for Observed Strata =====================================================\n');
            gui = plotFaciesFrequencyHistogram(gui, data);
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function CalculateAndPlotRandomModel_callback(source, eventsdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Calculate Random Models and Compare =====================================================\n');
            gui = calculateAndPlotRandomModels(gui, data); % 1 is plot flag, so plot the results
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function CalculateAndPlotOptimisedCycles_callback(source, eventsdata)
        
        if data.loaded
            fprintf('\n\n===================================================== Calculate Optimal Cycles =====================================================\n');
            gui = calculateAndPlotOptimisedCycles(gui, data); % 1 is plot flag, so plot the results
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function movingWindowAnalysisShowWindow_callback(source, eventsdata)
       if gui.movingWindowAnalysisWindowVisible == 0
            set(gui.movingWindow,'Visible','on');
            gui.movingWindowAnalysisWindowVisible = 1;
       end
    end

    function spectralAnalysisShowWindow_callback(source, eventsdata)
       if data.loaded
            fprintf('\n\n===================================================== Calculate Spectral Analysis =====================================================\n');
            gui = calculateSpectralAnalysis(gui, data); % 1 is plot flag, so plot the results
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No data loaded');
        end
    end

    function resetButton_callback(source, eventdata) 

        if isfield(gui,'f2') && ~isempty(gui.f2) && isvalid(gui.f2)
           close(gui.f2);
           gui = rmfield(gui,'f2');
        end

        if isfield(gui,'f3') && ~isempty(gui.f3) && isvalid(gui.f3)
           close(gui.f3);
           gui = rmfield(gui,'f3');
        end
        
        if isfield(gui,'f4') && ~isempty(gui.f4) && isvalid(gui.f4)
           close(gui.f4);
           gui = rmfield(gui,'f4');
        end

        if isfield(gui,'sp1') && ~isempty(gui.sp1) && isvalid(gui.sp1)
           cla(gui.sp1);
        end

        if isfield(gui,'sp2') && ~isempty(gui.sp2) && isvalid(gui.sp2)
           cla(gui.sp2);
        end

        if isfield(gui,'sp3') && ~isempty(gui.sp3) && isvalid(gui.sp3)
           cla(gui.sp3);
        end
     
        if isfield(gui,'movWindSectionPlot')
            cla(gui.movWindSectionPlot);
            close(gui.movWindSectionPlot);
            gui = rmfield(gui,'movWindSectionPlot');
        end
        
        if isfield(gui,'movWindPValuesPlot')
            cla(gui.movWindPValuesPlot);
            close(gui.movWindPValuesPlot);
            gui = rmfield(gui,'movWindPValuesPlot');
        end
        
        if gui.movingWindowAnalysisWindowVisible == 1
            set(gui.movingWindow,'Visible','off');
            gui.movingWindowAnalysisWindowVisible = 0;
        end

        data.loaded = 0;    % Boolean flag to show if data loaded or not; if not, subsequent functions disabled
        data.faciesCodes = zeros(1,1);
        data.faciesThick = zeros(1,1);
        data.sectLength = 0.0;
        data.maxNumbOfFacies = 0;
        data.sectionMean = 0.0;
    end

%% Now create the subordinate output window for the moving window analysis. Make invisible now, will be made visible when button in main GUI is cicked
    
    scrWidthProportion = 0.8;
    scrHeightIncrement = scrsz(4)/20; % Use this to space controls down right side of the main window
    controlStartY = (scrsz(4) * 0.80) - (scrHeightIncrement / 2);
    controlStartX = (scrsz(3) * scrWidthProportion) - 520;
    
    % position requires left bottom width height values. screensize vector
    % is in this format 1=left 2=bottom 3=width 4=height
    % Hide the window as it is constructed by setting visible to off
    gui.movingWindow = figure('Visible','off','Position',[1 1 scrsz(3)*scrWidthProportion scrsz(4)*0.8]);
    
   %  Construct the control panel components.
   hMovWindFPathLabel = uicontrol('style','text','string','MW analysis data file path:','Position',[controlStartX, controlStartY-scrHeightIncrement, 200, 15]);
   hMovWindFPath = uicontrol('Style','edit','String','sectionData/indusSuccessions/MWAnalysis/','Position',[controlStartX+200, controlStartY-scrHeightIncrement, 300, 25]);
%    hMovWindFPath = uicontrol('Style','edit','String','sectionData/','Position',[controlStartX+200, controlStartY-scrHeightIncrement, 200, 25]);
   
   hMovWindFNameLabel = uicontrol('style','text','string','MW analysis data filename:','Position',[controlStartX, (controlStartY-scrHeightIncrement*2), 200, 15]);
   hMovWindFName = uicontrol('Style','edit','String','utrankasSectionMWindowAnalysisBaseInc1SizeInc10.txt','Position',[controlStartX+200, (controlStartY-scrHeightIncrement*2), 300, 25]);
%     hMovWindFName = uicontrol('Style','edit','String','sectSwaps20recodedMWindowAnalysis.txt','Position',[controlStartX+200, (controlStartY-scrHeightIncrement*2), 200, 25]);

   hMovWindSizeLabel = uicontrol('style','text','string','MW analysis window size:','Position',[controlStartX, (controlStartY-scrHeightIncrement*3), 200, 15]);
   hMovWindSize = uicontrol('Style','edit','String','10','Position',[controlStartX+200, (controlStartY-scrHeightIncrement*3), 50, 25]);

   hMovWindPosLabel = uicontrol('style','text','string','MW analysis window base position:','Position',[controlStartX, (controlStartY-scrHeightIncrement*4), 200, 15]);
   hMovWindPos = uicontrol('Style','edit','String','1','Position',[controlStartX+200, (controlStartY-scrHeightIncrement*4), 50, 25]);

   hLoadMovWindData = uicontrol('Style','pushbutton','String','Load and plot moving window analysis data',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*5),200,25],...
          'BackgroundColor',[0.6 1.0 0.6],...
          'Callback',{@movingWindowAnalysisLoadAndPlot_callback});

   hPlotMovWindData = uicontrol('Style','pushbutton','String','Calculate and plot moving window analysis data',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*6),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@movingWindowAnalysisCalculateAndPlot_callback});
      
   hOneWindStats = uicontrol('Style','pushbutton','String','Calculate stats for one window',...
          'Position',[controlStartX+200,controlStartY-(scrHeightIncrement*7),200,25],...
          'BackgroundColor',[0.8 0.9 1.0],...
          'Callback',{@calculateAndPlotStatsForOneWindow_callback});
      
    function movingWindowAnalysisLoadAndPlot_callback(source, eventsdata)
        
        % This routine will produce a data output file, so need to get the file name stem for this and send it to the function
        dataFileNameAndPath = strcat(get(hMovWindFPath,'String'),get(hMovWindFName,'String'));
        
        if data.loaded
            fprintf('\n\n===================================================== Load & Plot Previous Moving Window Analysis =====================================================\n');
            gui = movingWindowAnalysisLoadAndPlot(gui, data, dataFileNameAndPath); 
        else
            message = sprintf('No data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No section data loaded');
        end
    end  
        
    function movingWindowAnalysisCalculateAndPlot_callback(source, eventsdata)
        
        % This routine will produce a data output file, so need to get the file name stem for this and send it to the function
        dataFileNameAndPath = strcat(get(hMovWindFPath,'String'),get(hMovWindFName,'String'));
        
        if data.loaded
            fprintf('\n\n===================================================== Calculate & Plot Moving Window Analysis =====================================================\n');
            gui = movingWindowAnalysisCalculateAndPlot(gui, data, dataFileNameAndPath); 
        else
            message = sprintf('No section data loaded\nUsed the load and plot section data button\n');
            m1 = msgbox(message,'No section data loaded');
        end
    end

    function calculateAndPlotStatsForOneWindow_callback(source, eventsdata)
        
         inputStr = get(hMovWindSize,'String');
         if isempty(str2num(inputStr))
            set(source,'string','0');
            warndlg('Input must be numerical');
         else
             windowSize = str2num(inputStr);
         end
         
         inputStr = get(hMovWindPos,'String');
         if isempty(str2num(inputStr))
            set(source,'string','0');
            warndlg('Input must be numerical');
         else
             windowPos = str2num(inputStr);
         end
        
         if data.loaded
             fprintf('\n\n===================================================== Calculate Stats for One Window =====================================================\n');
             gui = calculateAndPlotStatsForOneWindow(gui, data, windowSize, windowPos); 
         else
             message = sprintf('No section data loaded\nUsed the load and plot section data button\n');
             m1 = msgbox(message,'No section data loaded');
         end
    end
end

function [data, success] = loadSectionData(data, dataFileNameAndPath, coloursFileNameAndPath)
    
    % if we are loading a new section, need to clear out all the old
    % data, but cant just use clear because this gets rid of all arrays, so
    % need to think of a way to do this here...

    data.sectionName = dataFileNameAndPath;
    data.faciesCodes = zeros(1,1);
    data.thick = zeros(1,1);
    data.sectLength = 0.0;
    j = uint32(0);
    success = 1; % Assume loading will go ok, and reset to false below if not...

    % Read section data, thickness and facies
    if exist(dataFileNameAndPath, 'file')
        dataIn = load(dataFileNameAndPath);
        data.faciesThick = dataIn(:,1);
        data.faciesCodes = dataIn(:,2);
    else
        message = sprintf('Succession data file %s does not exist. Please correct the filename and try again.\n', dataFileNameAndPath);
        m1 = msgbox(message,'Data not loaded');
        success = 0;
        return;
    end
    
    checkMaxFaciesCode = max(data.faciesCodes);
    
    % Read colour map file, one colour and one text name per facies
    if exist(coloursFileNameAndPath, 'file');
        fid = fopen(coloursFileNameAndPath);
        checkCount = 0;
        while ~feof(fid) && checkCount <= checkMaxFaciesCode
            file = textscan(fid, '%u8 %f %f %f %s');
            checkCount = checkCount + 1;
        end
        fclose(fid);
        
    else
        message = sprintf('Facies colour map and name file %s does not exist. Please correct the filename and try again.\n', coloursFileNameAndPath);
        m1 = msgbox(message,'Data not loaded');
        success = 0;
        return;
    end
    
    % Extract facies names, numeric codes and RGB colours from data input file as a cell array
    data.faciesNames = file {5}(); % Because the 5th element in each row of the litho col file should be the facies name
    data.faciesNumbers = file {1}; % And the original faceis code is in 1st elemnt in each row - particularly important for optimal cycle analysis
    data.faciesColours = zeros(length(data.faciesNumbers),4); % Dimension colours matrix based on number of facies read from file
    for k=1:4
        data.faciesColours(:,k)= double(file{k}());
    end

    % Calculate the basic stats on the dataFacies array
    data.sectLength = max(size(data.faciesCodes));
    data.maxNumbOfFacies = size(data.faciesColours,1); % Return size of array dimension 1 which is number of facies. Dimension 2 will be 4 - code plus colours
    data.sectionMean = mean(data.faciesThick);

    % Make a reference section that uses the facies names to ID each unit
    % rather than the numerical facies code. Need to remember that faciesNames
    % is a cell array hence the curly brackets
    % Note that as of 12.7.2015 this array is not used elsewhere in the code
    for k=1:data.sectLength
        sectionFaciesNames{k} = data.faciesNames{data.faciesCodes(k)};
    end

    fprintf('Loaded vertical section data from %s\nFor %d units, total %d facies, mean unit thickness %4.3f m\n', dataFileNameAndPath, data.sectLength, data.maxNumbOfFacies, data.sectionMean);
end

function gui = plotObservedSectionInMainWindow(gui, data, bedCountMarker)

    % Plot the original data vertical section
    gui.sp1 = subplot('Position',[0.04 0.1 0.075 0.85]);
    cla(gui.sp1); % Clears the axis of whatever has been previously plotted
    
    [gui,~] = plotObservedSection(gui, data, 1, data.sectLength, bedCountMarker);
end
      
function [gui, vertSectAxes] = plotObservedSection(gui,data, sectionStart, sectionLength, bedCountMarker)
% general purpose routine to plot a section passed to the function in faciesCodes and faciesThick
% Note the full section data i
    hold on
    cumThick = 0;
    for j = sectionStart:sectionLength
        fCode = data.faciesCodes(j);
        yco = [cumThick, cumThick, cumThick + data.faciesThick(j), cumThick + data.faciesThick(j)];
        xco = [0, data.faciesCodes(j), data.faciesCodes(j), 0];
        faciesCol = [data.faciesColours(fCode,2) data.faciesColours(fCode,3) data.faciesColours(fCode,4)];
        patch(xco, yco, faciesCol,'EdgeColor','none');
        cumThick = cumThick + data.faciesThick(j);
    end
    
    plotBedCountMarkers(data, bedCountMarker);
    
    grid on;
    set(gca,'Layer','top');
    xlabel('Facies code');
    ylabel('Thickness (m)');
    set(gca,'XTick', 1:(data.maxNumbOfFacies),'XTickLabel', data.faciesNames, 'TickDir', 'out', 'XTickLabelRotation',90);
    hold off;
    vertSectAxes = gca;
end

function plotBedCountMarkers(data, bedCountMarker)

    % Plot bed count markers, 1 line every bedCountMarker number of beds
    for j = 1:bedCountMarker:data.sectLength
        xco = [data.maxNumbOfFacies, data.maxNumbOfFacies + 1];
        totalThickSoFar = sum(data.faciesThick(1:j));
        yco = [totalThickSoFar, totalThickSoFar];
        line(xco, yco);
    end
end

function gui = plotObservedSectionInScrollWindow(gui, data, bedCountMarker)

    % Plot the original data vertical section
    scrsz = get(0,'ScreenSize'); % screen dimensions vector
    gui.sectScrollFig = figure('Visible','on','Position',[1 ,1, scrsz(3)/4, scrsz(4)]);
    [gui, vertSectAxes] = plotObservedSection(gui, data, 1, data.sectLength, bedCountMarker);
    
    plotBedCountMarkers(data, bedCountMarker);
    
    % Set axis limits 
    maxY = sum(data.faciesThick);
    dy = 5; % Specify scroll increment
    set(vertSectAxes,'xlim',[0, data.maxNumbOfFacies + 1]); % facies number plus room for bed count lines
    set(vertSectAxes,'ylim',[0, maxY]);
    
    % define slider position
    pos = get(vertSectAxes,'position');
    newPos = [pos(3) + 0.15, pos(2), pos(3) + 0.15, pos(4)]; % define slider position just right of axis
    
    % Setting up callback string to modify yLim of axis (gca) based on the position of the slider (gcbo)
    sliderCallBack = ['set(gca,''ylim'',get(gcbo,''value'')+[0 ' num2str(dy) '])'];
    uicontrol('style','slider','units','normalized', 'position', newPos, 'callback',sliderCallBack, 'min',0, 'max',maxY);
    set(gcf,'doublebuffer','on'); % This avoids flickering when updating the axis
   
end

function gui = plotBedThicknessBarChart(gui, data, bedCountMarker)

    % Plot the bed thickness as a bar chart
    scrsz = get(0,'ScreenSize'); % screen dimensions vector
    gui.bedThicknessBarChartFig = figure('Visible','on','Position',[1 ,1, scrsz(3)*0.75, scrsz(4)/2]);
    h1 = bar(data.faciesThick,'FaceColor','flat'); % Flat face colour so it can be changed in the next command
%     for j = 1:bedCountMarker:numel(data.faciesThick)
%         h1.CData(j,:) = [1,0,0]; % Set every bedCountMarker th bed to red colour
%     end
    
    h1.CData(:,:) = [data.faciesColours(data.faciesCodes,2), data.faciesColours(data.faciesCodes,3), data.faciesColours(data.faciesCodes,4)]; % Set each bar to corresponding facies colour
    xlabel("Layer number");
    ylabel("Layer thickness (m)");
    grid on;
    
      % Plot the bed thickness as a bar chart with a log thickness scale
    scrsz = get(0,'ScreenSize'); % screen dimensions vector
    gui.bedThicknessBarChartFig = figure('Visible','on','Position',[50 ,50, scrsz(3)*0.75, scrsz(4)/2]);
    logBedThicknessData = log10(data.faciesThick);
    h2 = bar(logBedThicknessData,'FaceColor','flat');
%     for j = 1:bedCountMarker:numel(data.faciesThick)
%         h2.CData(j,:) = [1,0,0]; % Set every bedCountMarker th bed to red colour
%     end
    h2.CData(:,:) = [data.faciesColours(data.faciesCodes,2), data.faciesColours(data.faciesCodes,3), data.faciesColours(data.faciesCodes,4)];
    xlabel("Layer number");
    ylabel("Log layer thickness (m)");
    grid on;
    
end

function [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(faciesCodes, faciesNames, plotFlag, g1)

    % find the number of elements in the succession. NB size returns both matrix dimensions so max ensures nz is the biggest which should be the length of the section
    nz = max(size(faciesCodes));
    
    % Find the maximum facies code used in the facies succession - this is the size for both dimensions of the TP matrix which can now be defined
    nFacies = max(faciesCodes);
    
    m = 0; % number of transitions, not really needed but never mind
    TFMatrix = zeros(nFacies, nFacies);
    TPMatrix = zeros(nFacies, nFacies);
    TPMultiplierMatrix = zeros(nFacies, nFacies);
    
    % populate the multiplier matrix with value 10 on the 1-offset diagonal - where the max p values will be for ABCDE strata
    for i = 1:nFacies-1
        TPMultiplierMatrix(i,i+1) = 10.0;
    end
    TPMultiplierMatrix(nFacies,1) = 10.0; % Complete the diagonal with the "wraparound" position in the matrix corner
    
    % Now loop through the elements in the succession and for each different facies from-to transition, increment the appropriate cell in the matrix
    for i =1 : nz-1
        fromFacies = faciesCodes(i);
        toFacies = faciesCodes(i+1);
        % mark transitions between different facies
        if fromFacies > 0 && toFacies > 0 && fromFacies ~= toFacies % Make sure facies codes are not zero because zero values would record an error
            TFMatrix(fromFacies, toFacies) = TFMatrix(fromFacies, toFacies) + 1; % increment the appropriate value in the tp matrix
            m = m + 1;
        end     
    end

    % Now calculate the transition probability matrix from the transition frequency matrix
    rowSums=sum(TFMatrix,2); % Calculates the sum of each row in TF matrix and stores as vector rowSums
    for k=1:nFacies
        for j=1:nFacies
            if rowSums(k) > 0 % if rowsum > 0 divide TF value by row sum to get transition probability
                TPMatrix(k,j) = TFMatrix(k,j) / rowSums(k);
            else
                TPMatrix(k,j) = 0;
            end
        end
    end
    
    % Now calculate the Markov order metrics

    % Calculate the metric for the maximum diagonal
    diagMetric = zeros(1,nFacies-1);
    diagMetricMultiplied = zeros(1,nFacies-1);
    TPMatrixMultiplied = TPMultiplierMatrix .* TPMatrix;
    
    % Now loop through each offset diagonal in the matrix and calculate the
    % mean of the cell values in the diagonal
    % Also do this for the version of the TP matrix with a *10 weighted one-offset diagonal
    for j=1:nFacies-1
    
        diagMetric(j) = (sum(diag(TPMatrix,j)) + sum(diag(TPMatrix,-(nFacies-j))) )/ nFacies; % calculate a mean for the jth FULL diagonal in the matrix
        diagMetricMultiplied(j) = (sum(diag(TPMatrixMultiplied,j)) + sum(diag(TPMatrixMultiplied,-(nFacies-j))) )/ nFacies; % Also calculate a mean based on the weighted one-offset diagonal
    end
    
    markovOrderMetric = max(diagMetric)- min(diagMetric);
    oneOffsetMValueDiagCheck = max(diagMetricMultiplied)- min(diagMetricMultiplied);

    if plotFlag == 1
        
        TPCellSize = 1.0;
        matrixTopYco = (nFacies + 1) * TPCellSize;
        matrixBottYco = 0.0;
        
        for i = 1:nFacies
            for j = 1:nFacies
                yco = [matrixBottYco+((j-0.5)*TPCellSize) matrixBottYco+((j-0.5)*TPCellSize) matrixBottYco+((j+0.5)*TPCellSize) matrixBottYco+((j+0.5)*TPCellSize)];
                xco = [(i-0.5)*TPCellSize (i+0.5)*TPCellSize (i+0.5)*TPCellSize (i-0.5)*TPCellSize];
                labelStr = sprintf('%3.2f',TPMatrix(j,i));

                if TPMatrix(j,i) <= 0.5
                    colVect = [1 TPMatrix(j,i)/0.5 TPMatrix(j,i)/1.1111111]; % gradation of colours from red (p=0) to orange-yellow (p=0.5)
                else
                    colVect = [1-(TPMatrix(j,i)) 1 1-(TPMatrix(j,i)/1.111111)]; % gradation of colours from orange-yellow (p=0.5) to green (p=1)
                end
               
                if i==j
                    colVect = [0.9 0.9 0.9]; % Matrix diagonal is a special case so colour light grey
                end
                
                patch(xco, yco, colVect);
                
                if i ~= j % So do not add text labels to the diagonal
                  text(double(((i-0.7)*TPCellSize))+(TPCellSize*0.33), double(matrixBottYco+(j*TPCellSize)-0.3), labelStr,'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
                end
            end
        end
        
        labelStr = cell(1,nFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
        for k=1:nFacies
            labelStr{k} = faciesNames{k}; % Because bestFaciesOrder is not in the same order as faciesNames
        end
        g1.XTickLabel = labelStr;
        g1.YTickLabel = labelStr;
        
        set(g1,'XTick',1:nFacies);
        set(g1,'YTick',1:nFacies);
       
        xlabel('Facies Code: To', 'FontSize',10);
        ylabel('Facies code: From', 'FontSize',10);
        axis tight;
    end
end

function gui = plotFaciesFrequencyHistogram(gui, data)
% Plot the histogram of the facies frequencies

    gui.sp3 = subplot('Position',[0.2 0.075 0.3 0.4]);
    cla(gui.sp3);
    hold on
    
    frequencyDataFacies=histc(data.faciesCodes,1:data.maxNumbOfFacies);
    for j=1:data.maxNumbOfFacies
        xco = [j-1, j-1, j, j];
        yco = [0, frequencyDataFacies(j), frequencyDataFacies(j), 0];
        faciesCol = [data.faciesColours(j,2) data.faciesColours(j,3) data.faciesColours(j,4)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
    end

    set(gca,'XTick', 0.5:(data.maxNumbOfFacies-0.5),'XTickLabel', data.faciesNames);
    xlabel('Lithofacies');
    ylabel('Frequency');
    
    hold off;
    
    % Check is all facies are present in the section, and if not, issue a warning to recode
    zeroFreqWarningFlag = 0;
    for i=1:data.maxNumbOfFacies
        if frequencyDataFacies(i) == 0
             zeroFreqWarningFlag = 1;
        end
    end
    if zeroFreqWarningFlag
        message = sprintf('At least one facies code has no occurences in this section\nConsider re-coding facies with concurruent numbering to avoid misleading results.');
        m1 = msgbox(message,'Missing Facies');
    end 
end

function gui = calculateAndPlotRandomModels(gui, data)

    TRUE = uint8(1);
    FALSE = uint8(0);
    maxIterations = 5000;
    numberOfSwaps = data.sectLength;
    maxRun = 3;
    minRun = 0;
    runBinIncrement = 0.05;
    runRange = maxRun - minRun;

    % Calculate and output the order metric for the entire data succession
    [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    runsOrderMetric = calculateRunsOrderMetric(data.faciesThick);
    fprintf('Markov metric for strata is %4.3f\nRuns analysis metric for strata is %5.4f\n', markovOrderMetric, runsOrderMetric);

    % Now calculate the metrics for many iterations of a random model

    % Shuffle the observed section and calculate the facies and thickness order metrics each time
    for j = 1:maxIterations;
        if data.maxNumbOfFacies > 3
            [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(data.faciesCodes, data.faciesThick, numberOfSwaps);
            multiMarkovOrderMetricDataShuffled(j) = calculateTPMatrixAndOrderMetric(shuffledFacies, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01])); % note the zero is the TP matrix plot flag, set to false so no plot
        else
            shuffledThick = shuffleSectionThicknessOnly(data.faciesThick, numberOfSwaps);
        end
        
        multiRunsOrderMetricDataShuffled(j) = calculateRunsOrderMetric(shuffledThick); 
    end
    
    fprintf('For %d iterations of a SHUFFLED DATA model\n', maxIterations);

    % Stats on the shuffled section random model - facies
    if data.maxNumbOfFacies > 3
        meanMultiMarkovDataShuffled = mean(multiMarkovOrderMetricDataShuffled);
        stdDevMultiMarkovDataShuffled = std(multiMarkovOrderMetricDataShuffled);
        bins = 0:0.02:1.00; % because 0<m<=1
        multiMarkovOrderMetricDataShuffledHistData=histc(multiMarkovOrderMetricDataShuffled, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
        markovPValueSum = sum(multiMarkovOrderMetricDataShuffledHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1
        fprintf('Markov stats mean %5.4f std dev %5.4f Markov order metric P Value %5.4f\n', meanMultiMarkovDataShuffled, stdDevMultiMarkovDataShuffled, markovPValueSum);
    else
        fprintf('Only %d distinct facies, minimum of 4 needed to shuffle sections and run Markov stats\n', data.maxNumbOfFacies);
    end
        
    % Stats on the shuffled section random model - thickness
    meanMultiRunsDataShuffled = mean(multiRunsOrderMetricDataShuffled);
    stdDevMultiRunsDataShuffled = std(multiRunsOrderMetricDataShuffled);
    bins = 0:0.05:runRange; % 0 is the minimum run metric, 3 a generally maximum value (this is what runRange should be set to)
    multiRunsOrderMetricDataShuffledHistData = histc(multiRunsOrderMetricDataShuffled, bins)/maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    runBinIndex = 1 + int16(((runsOrderMetric-minRun)/runRange)*(runRange/runBinIncrement)); % Position of runs stat in the histogram
    runsPValueSum = sum(multiRunsOrderMetricDataShuffledHistData(runBinIndex:length(multiRunsOrderMetricDataShuffledHistData))); % area under curve from r to max run value
    fprintf('Runs stats mean %5.4f std dev %5.4f Runs analysis metric P Value %5.4f\n', meanMultiRunsDataShuffled, stdDevMultiRunsDataShuffled, runsPValueSum);

    % Plot the results
    scrsz = get(0,'ScreenSize'); % screen dimensions vector

    % Plot the original vertical section, alongside the facies and thickness
    % order plots showing the outcrop datapoint and the multiple iteration shuffled section
    % metric frequency distribution
    gui.f2 = figure('Visible','on','Position',[1 scrsz(4)/4 (scrsz(3)/2) (scrsz(4)/3)*2]);
    
    subplot('Position',[0.075 0.1 0.075 0.85]);
    plotShuffledSection(data);
    
    if data.maxNumbOfFacies > 3
        plotMCMarkovHistogram(multiMarkovOrderMetricDataShuffledHistData, markovOrderMetric, [0.25 0.57 0.65 0.38]);
    end
    
    plotMCRunsHistogram(minRun, runBinIncrement, maxRun, multiRunsOrderMetricDataShuffledHistData, runsOrderMetric);
    
    % Complete analysis with a textbox message summary
    if data.maxNumbOfFacies > 3
        message = sprintf('Markov metric %4.3f\nRuns analysis metric %5.4f\nMarkov order metric P Value %5.4f\nRuns analysis metric P Value %5.4f\n', markovOrderMetric, runsOrderMetric, markovPValueSum, runsPValueSum);
    else
        message = sprintf('Runs analysis metric %5.4f\nRuns analysis metric P Value %5.4f\n', runsOrderMetric, runsPValueSum);
    end
    m1 = msgbox(message,'Random model comparison');
end

function plotMCMarkovAndRunsHistograms(data, markovOrderMetric, multiMarkovOrderMetricDataShuffledHistData, runsOrderMetric, minRun, runBinIncrement, maxRun, multiRunsOrderMetricDataShuffledHistData, markovPValueSum, runsPValueSum)

    plotMCMarkovHistogram(multiMarkovOrderMetricDataShuffledHistData, markovOrderMetric, [0.25 0.57 0.65 0.38]);

    plotMCRunsHistogram(minRun, runBinIncrement, maxRun, multiRunsOrderMetricDataShuffledHistData, runsOrderMetric);
    
    % Complete analysis with a textbox message summary
    message = sprintf('Markov metric %4.3f\nRuns analysis metric %5.4f\nMarkov order metric P Value %5.4f\nRuns analysis metric P Value %5.4f\n', markovOrderMetric, runsOrderMetric, markovPValueSum, runsPValueSum);
    m1 = msgbox(message,'Random model comparison');
end

function  plotMCMarkovHistogram(multiMarkovOrderMetricDataShuffledHistData, markovOrderMetric, plotCoords)
    % Subplot for the Markov order analysis histogram   
    sp1 = subplot('Position', plotCoords);
    cla(sp1); % Clears the axis of whatever has been previously plotted
    hold on
    bins = 0:0.02:1.00; % Make sure bins is set correctly for Markov plots
    bar(bins, multiMarkovOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]); % Colour is dark slate blue
    maxFreq = max(multiMarkovOrderMetricDataShuffledHistData) *1.1; % This is needed to scale the plot
    x = [markovOrderMetric markovOrderMetric];
    y = [0 maxFreq];
    line(x,y, 'color', [0.80 0.00 0.00], 'linewidth', 3.0); % Colour is dark red
    grid on;
    axis([0 1 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Markov Order Metric for Facies');
    ylabel('Relative Frequency');
end

 function plotMCRunsHistogram(minRun, runBinIncrement, maxRun, multiRunsOrderMetricDataShuffledHistData, runsOrderMetric)
    % Subplot for the runs analysis histogram    
    subplot('Position',[0.25 0.1 0.65 0.38]);
    hold on
    bins = minRun:runBinIncrement:maxRun; % Make sure that bins is set correctly for Runs analysis plots
    bar(bins, multiRunsOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]);
    maxFreq = max(multiRunsOrderMetricDataShuffledHistData) * 1.1; % This is needed to scale the plot

    lineColor = [0.80 0.00 0.00];
    x = [runsOrderMetric runsOrderMetric];
    y = [0 maxFreq]; % Draw the data line from y=0 to y=max frequency of the three histograms
    line(x,y, 'color', lineColor, 'linewidth', 3.0);

    grid on;
    axis([0 Inf 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Runs Analysis Order Metric for Thickness ');
    ylabel('Relative Frequency');
end

% function plotShuffledSection(gui, data)
function plotShuffledSection(data)

    % Plot a shuffled iteration of the observed vertical section
    numberOfSwaps = data.sectLength;
    [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(data.faciesCodes, data.faciesThick, numberOfSwaps);
    cumThick = 0;
    for j=1:data.sectLength
        fCode = shuffledFacies(j);
        yco = [cumThick, cumThick, cumThick + shuffledThick(j), cumThick + shuffledThick(j)];
        xco = [0, shuffledFacies(j), shuffledFacies(j), 0];
        faciesCol = [data.faciesColours(fCode,2) data.faciesColours(fCode,3) data.faciesColours(fCode,4)];
        patch(xco, yco, faciesCol,'EdgeColor','none');
        cumThick = cumThick + shuffledThick(j);
    end
    grid on;
    set(gca,'Layer','top');
    xlabel('Facies code');
    ylabel('Thickness (m)');
    set(gca,'XTick', 1:(data.maxNumbOfFacies),'XTickLabel', data.faciesNames, 'TickDir', 'out', 'XTickLabelRotation',90);
    title('Shuffled section');
end

function [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(sectFacies, sectThick, totalSwaps)
% function to shuffle the facies and thickness succession to ensure a random configuration, and assuring no adjacent same facies occurrences in final section
% NB only works if more than three distinct facies in the section because not possible to avoid adjacent facies occurrences with fewer

    % Make copies of the original data in new arrays that will be used to
    % store the shuffled sections
    shuffledFacies = sectFacies;
    shuffledThick = sectThick;
    n = uint16(max(size(shuffledFacies)));
    j = 0;
    infiniteStopCount = 0;
    
    while j < totalSwaps && infiniteStopCount < 1000000
        
        % Select two unit numbers randomly to be swapped
        unit1 = uint16((rand * (n-1)) + 1);
        unit2 = uint16((rand * (n-1)) + 1);
        
        % Need to check above and below for both positions that swapping will not put same
        % facies adjacent to one another and cause a transition to self
        swapFacies1 = shuffledFacies(unit1);
        if unit1 > 1 swapFacies1Below = shuffledFacies(unit1-1); else swapFacies1Below = 0;end
        if unit1 < n swapFacies1Above = shuffledFacies(unit1+1); else swapFacies1Above = 0;end
        
        swapFacies2 = shuffledFacies(unit2);
        if unit2 > 1 swapFacies2Below = shuffledFacies(unit2-1); else swapFacies2Below = 0;end
        if unit2 < n swapFacies2Above = shuffledFacies(unit2+1); else swapFacies2Above = 0;end
        
        % So compare facies in their new positions with the facies above and below and
        % only swap and increment loop counter if NOT the same...
        if swapFacies1Below ~= swapFacies2 && swapFacies1Above ~= swapFacies2 && swapFacies2Below ~= swapFacies1 && swapFacies2Above ~= swapFacies1
            
            %Swap the facies
            temp = shuffledFacies(unit1);
            shuffledFacies(unit1) = shuffledFacies(unit2);
            shuffledFacies(unit2) = temp;

            %Swap the thicknesses
            temp = shuffledThick(unit1);
            shuffledThick(unit1) = shuffledThick(unit2);
            shuffledThick(unit2) = temp;

            j = j + 1;
        end
        
        infiniteStopCount = infiniteStopCount + 1;
    end
end

function shuffledThick = shuffleSectionThicknessOnly(sectThick, totalSwaps)
% function to shuffle the thickness succession only, to ensure a random configuration

    % Make copies of the original data in new arrays that will be used to store the shuffled section
    shuffledThick = sectThick;
    n = uint16(max(size(shuffledThick)));
    j = 0;
    
    while j < totalSwaps
        
        % Select two unit unit numbers randomly to be swapped
        unit1 = uint16((rand * (n-1)) + 1);
        unit2 = uint16((rand * (n-1)) + 1);
        
        %Swap the thicknesses
        temp = shuffledThick(unit1);
        shuffledThick(unit1) = shuffledThick(unit2);
        shuffledThick(unit2) = temp;

        j = j + 1;
    end
end

function runsOrderMetric = calculateRunsOrderMetric(thicknesses)

    % find the number of units in the succession and declare arrays accordingly
    nz = max(size(thicknesses));
    deltaThick = zeros(1,nz);
    runsUp = zeros(1,nz);
    runsDown = zeros(1,nz);

    % Calculate the change in thickness between successive units
    i = 1:nz-1;
    j =2:nz; % so j = i + 1 therefore thickness change is thickness(j) - thickness(i)
    deltaThick(i) = thicknesses(j) - thicknesses(i);

    if deltaThick(1) > 0 runsUp(1) = 1; end
    if deltaThick(1) < 0 runsDown(1) = 1; end

    for i=2:nz
        if deltaThick(i) > 0 runsUp(i) = runsUp(i-1)+1; end
        if deltaThick(i) < 0 runsDown(i) = runsDown(i-1)+1; end
    end

    runsUpNormSum = (sum(runsUp)/nz);
    runsDownNormSum = (sum(runsDown)/nz);
    runsOrderMetric = (runsUpNormSum + runsDownNormSum);
end

function gui = calculateAndPlotOptimisedCycles(gui, data) % 1 is plot flag, so plot the results

    TRUE = uint8(1);
    FALSE = uint8(0);
    ALL_PLOTS = TRUE; % TRUE if you want a TP matrix and ideal cycle plot for each of the optimum m and diagonal permutations. FALSE if you just want summary plots
    results.testFaciesSect = 0; % results is passed to the functions below and returned full of results, but needs to be created here to be passed
    
    tic

    % Calculate the basic stats on the dataFacies array
    sectionLength = max(size(data.faciesCodes));
    maxNumbOfFacies = length(data.faciesColours);
    sectionMean = mean(data.faciesThick);

    fprintf('For %d units, total %d facies, mean unit thickness %4.3f m\n', sectionLength, maxNumbOfFacies, sectionMean);

    % Calculate and output the order metric for the entire data succession
    markovOrderMetric = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    fprintf('Markov metric for strata is %5.4f\n', markovOrderMetric);

    
    
    % Calculate the various statistics returned from the function for each facies code permutation
    results = calcStatsForEachFaciesOrderPermutation(data, sectionLength, results); 
   
    % Now loop again to find and record all of the permutations and resulting strat sections that gave the maximum or minimum m and diag values found above
    % and record various information about each permutation
    results = findOptimumAndWorstPermutations(data, sectionLength, maxNumbOfFacies, results);
    
    toc % Gives time taken to do all the permutation calculations

    % **************************************************************************************
    % Report all of the results as text output
    fullFileName = strcat('optimalCyles', '_report.txt'); % Save the permutation results to .mat file in ascii format
    fOut = fopen(fullFileName,'w');
    fprintf('Largest markov stat value %5.4f found from %d permutations \n', results.maxMarkov, results.maxMarkovCount);
    fprintf(fOut,'Largest markov stat value %5.4f found from %d permutations \n', results.maxMarkov, results.maxMarkovCount);
    fprintf('Lowest markov stat value %5.4f found from %d permutations \n', results.minMarkov, results.minMarkovCount);
    fprintf(fOut,'Lowest markov stat value %5.4f found from %d permutations \n', results.minMarkov, results.minMarkovCount);
    fprintf('Highest m value * one-offset diagonal value %5.4f found from %d permutations \n', results.maxMarkovDiagProduct, results.maxDiagCount);
    fprintf(fOut,'Highest m value * one-offset diagonal value %5.4f found from %d permutations \n', results.maxMarkovDiagProduct, results.maxDiagCount);
    fprintf('Optimisation indicator: %5.4f From %d total permutations, %d permutations gave optimal score\n', results.maxDiagCount/results.permsCount, results.permsCount, results.maxDiagCount);
    fprintf(fOut,'Optimisation indicator: %5.4f From %d total permutations, %d permutations gave optimal score\n', results.maxDiagCount/results.permsCount, results.permsCount, results.maxDiagCount);
    
    k=1:maxNumbOfFacies; % NB implicit loop on k
    fprintf('Original facies: ');
    fprintf('%s ', data.faciesNames{k});
    fprintf('\nOptimal codings are:');
    fprintf(fOut,'Original facies: ');
    fprintf(fOut,'%s ', data.faciesNames{k});
    fprintf(fOut,'\nOptimal codings are:');
    for j=1:results.maxDiagCount
         fprintf('\n%d ',j);
         fprintf(fOut,'\n%d ',j); 
         fprintf('%s ', data.faciesNames{results.bestDiagFaciesOrder(j,k)}); % NB the implied loop on k in this - separate formatting fprints required above therefore
         fprintf(fOut,'%s ',data.faciesNames{results.bestDiagFaciesOrder(j,k)});
         fprintf('m=%5.4f diag=%5.4f m-diag product %5.4f', results.markovOrderMetricTest(results.maxDiagPerm(j)), results.oneOffsetDiagCheck(results.maxDiagPerm(j)), results.markovOrderMetricTest(results.maxDiagPerm(j)) * results.oneOffsetDiagCheck(results.maxDiagPerm(j)) );
         fprintf(fOut,'m=%5.4f diag=%5.4f m-diag product %5.4f', results.markovOrderMetricTest(results.maxDiagPerm(j)), results.oneOffsetDiagCheck(results.maxDiagPerm(j)), results.markovOrderMetricTest(results.maxDiagPerm(j)) * results.oneOffsetDiagCheck(results.maxDiagPerm(j)) );
     end

    fprintf('\n');
    fprintf(fOut,'\n');
    
    k=1:maxNumbOfFacies; % NB implicit loop on k
    fprintf('\nWorst codings are:');
    fprintf(fOut,'\nWorst codings are:');
    for j=1:results.minMarkovCount
        fprintf('\n%d ',j);
        fprintf(fOut,'\n%d ',j);
        fprintf('%s ', data.faciesNames{results.worstFaciesOrder(j,k)}); % NB the implied loop on k in this - separate formatting fprints required above therefore
        fprintf(fOut,'%s ', data.faciesNames{results.worstFaciesOrder(j,k)});
    end
    fprintf('\n');
    fprintf(fOut,'\n');
    fclose(fOut);

    fullFileName = strcat('optimalCycles', '_report_perms.txt'); % Save the permutation results to .mat file in ascii format
    dataDump = results.markovOrderMetricTest;
    save(fullFileName, 'dataDump','-ascii');
      
    % Now finally, compare the first optimally coded section with random
    % shuffled sections to see if the m value gives any evidence of an
    % orderd succession
    maxIterations = 5000;
    numberOfSwaps = data.sectLength;

    % Calculate and output the order metric for the first best diagonal optimally coded succession
    [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(results.bestDiagFaciesSect(1,:), data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));  % note the zero is the TP matrix plot flag, set to false to not trigger a plot
    fprintf('Markov metric for optimised facies order strata is %4.3f\n', markovOrderMetric);

    % Now calculate the metrics for many iterations of a random model
    % Shuffle the observed section and calculate the order metric and DW stat each time
    for j = 1:maxIterations;
        [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(data.faciesCodes, data.faciesThick, numberOfSwaps);
        multiMarkovOrderMetricDataShuffled(j) = calculateTPMatrixAndOrderMetric(shuffledFacies, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01])); % note the zero is the TP matrix plot flag, set to false so no plot
    end

    % Stats on the shuffled section random model
    meanMultiMarkovDataShuffled = mean(multiMarkovOrderMetricDataShuffled);
    stdDevMultiMarkovDataShuffled = std(multiMarkovOrderMetricDataShuffled);
    bins = 0:0.02:1.00; % because 0<m<=1
    multiMarkovOrderMetricDataShuffledHistData=histc(multiMarkovOrderMetricDataShuffled, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
    markovPValueSum = sum(multiMarkovOrderMetricDataShuffledHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1

    fprintf('For %d iterations of a SHUFFLED DATA model\n', maxIterations);
    fprintf('Markov stats mean %5.4f std dev %5.4f For optimised facies coding Markov order metric P Value %5.4f\n', meanMultiMarkovDataShuffled, stdDevMultiMarkovDataShuffled, markovPValueSum);
    
% *************************************************************************************
    % Plot the results graphically
    % scrsz = left bottom width height
    scrsz = get(0,'ScreenSize'); % screen dimensions vector
    gui.f3 = figure('Visible','on','Position',[1 1 (scrsz(3)*0.75) (scrsz(4)*0.75)]);

    % TP matrix for the original facies coding
    h1 = subplot('Position', [0.05 0.55 0.20 0.3]);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(data.faciesCodes, data.faciesNames, 1, h1);
    axis tight;
    labelStr = cell(maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    h1.XTickLabel = {data.faciesNames{1:maxNumbOfFacies}}; % because this is a plot for the original order of facies row coding, just use that
    h1.YTickLabel = {data.faciesNames{1:maxNumbOfFacies}};
    str = sprintf('TP Matrix, original coding\nm=%4.3f diag=%4.3f m*diag %4.3f', mStat, diag, mStat * diag);
    title(str);

    % Plot the original facies coding order
    h2 = subplot('Position',[0.27 0.55 0.03 0.3]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [data.faciesColours(k,2) data.faciesColours(k,3) data.faciesColours(k,4)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
        labelStr = data.faciesNames(k);
        text(1.25, k+0.2, labelStr, 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end
    
    % TP matrix for the worst permutation
    h3 = subplot('Position',[0.35 0.55 0.20 0.3]);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(results.worstFaciesSect(1,:), data.faciesNames, 1, h3);
    axis tight;
    colourCode = zeros(maxNumbOfFacies, 3); % make sure colour code matrix is ready but empty
    labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    for k=1:maxNumbOfFacies
        % For each jth permutation in each worstFaciesOrder the array elements are in original facies order, 
        % so element j,k gives the facies k row code for the jth worst permutation, in this case worst permutation 1
        labelStr{results.worstFaciesOrder(1,k)} = data.faciesNames{k}; 
        colourCode(results.worstFaciesOrder(1,k),1:3) = data.faciesColours(k,2:4);
    end
    h3.XTickLabel = labelStr;
    h3.YTickLabel = labelStr;
    str = sprintf('Worst TP Matrix, worst facies coding\nm=%4.3f m*diags %4.3f', mStat, mStat * diag);
    title(str);
    
    % Plot a worst permutation facies coding order
    h4 = subplot('Position',[0.57 0.55 0.03 0.3]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)]; % Colour codes are set in the TP plot code above
        patch(xco, yco, faciesCol,'EdgeColor','black');
        text(1.25, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end
    
    % TP matrix for the first best permutation
    h5 = subplot('Position',[0.65 0.55 0.20 0.3]);
    [mStat,diag] = calculateTPMatrixAndOrderMetric(results.bestDiagFaciesSect(1,:), data.faciesNames, 1, h5);
    axis tight;
    colourCode = zeros(maxNumbOfFacies, 3); % make sure colour code matrix is ready but empty
    labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
    for k=1:maxNumbOfFacies
        % For each jth permutation in each bestDiagFaciesOrder the array elements are in original facies order, 
        % so element j,k gives the facies k row code for the jth best permutation, in this case best permutation 1
        labelStr{results.bestDiagFaciesOrder(1,k)} = data.faciesNames{k}; 
        colourCode(results.bestDiagFaciesOrder(1,k),1:3) = data.faciesColours(k,2:4);
    end
    h5.XTickLabel = labelStr;
    h5.YTickLabel = labelStr;
    str = sprintf('Best TP Matrix, optimal facies coding\nm=%4.3f m*diags %4.3f', mStat, mStat * diag);
    title(str);

    % Plot a best permutation facies coding order
    h6 = subplot('Position',[0.87 0.55 0.03 0.3]);
    for k=1:maxNumbOfFacies
        yco = [k k k+1 k+1];
        xco = [1 2 2 1];
        faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)];
        patch(xco, yco, faciesCol,'EdgeColor','black');
        text(1.25, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
        axis off;
    end

     % plot the histogram of m values from all permutations
    h4 = subplot('Position',[0.05 0.1 0.18 0.35]);
    histogram(results.markovOrderMetricTest, 'EdgeColor','none', 'FaceColor','blue');
    xlabel('m value', 'FontSize',10);
    ylabel('Frequency', 'FontSize',10);
    grid on;
    
    % plot the histogram of diagonal scores from the n maximum score permutations
    h5 = subplot('Position',[0.27 0.1 0.18 0.35]);
    histogram(results.markovDiagProduct, 'EdgeColor','none', 'FaceColor','blue');
    xlabel('d', 'FontSize',10);
    ylabel('Frequency', 'FontSize',10);
    grid on;

    % Subplot for the Markov order analysis histogram   
    h6 = subplot('Position',[0.50 0.1 0.40 0.35]);
    bins = 0:0.02:1.00; % Make sure bins is set correctly for Markov plots
    bar(bins, multiMarkovOrderMetricDataShuffledHistData, 'EdgeColor','none', 'BarWidth', 1, 'FaceColor',[0.2 0.4 0.7]); % Colour is dark slate blue
    maxFreq = max(multiMarkovOrderMetricDataShuffledHistData) *1.1; % This is needed to scale the plot
    x = [markovOrderMetric markovOrderMetric];
    y = [0 maxFreq];
    line(x,y, 'color', [0.80 0.00 0.00], 'linewidth', 3.0); % Colour is dark red
    grid on;
    axis([0 1 0 Inf]);
    set(gca,'Layer','top');
    xlabel('Markov Order Metric for Facies');
    ylabel('Relative Frequency');
    
    % If flag set, plot all of the best M-value permutation TP matrices and corresponding facies order summary sections
    if ALL_PLOTS == TRUE
        
        scrsz = get(0,'ScreenSize'); % screen dimensions vector
     
        for j = 1:results.maxMarkovCount
            
            figure('Visible','on','Position',[scrsz(3)/j (scrsz(4)*0.5) (scrsz(3)*0.2) (scrsz(4)*0.2)]);
            hX = subplot('Position',[0.05 0.10 0.70 0.8]);
            % TP matrix for the jth best M-value permutation
            [mStat,diag] = calculateTPMatrixAndOrderMetric(results.bestFaciesSect(j,:), data.faciesNames, 1, h5);
            axis tight;
            colourCode = zeros(maxNumbOfFacies, 3); % make sure colour code matrix is ready but empty
            labelStr = cell(1,maxNumbOfFacies); % Labels need to be in a cell array in order to put in X and YTickLabel
            for k=1:maxNumbOfFacies
                % For each jth permutation in each bestDiagFaciesOrder the array elements are in original facies order, 
                % so element j,k gives the facies k row code for the jth best permutation, in this case best permutation 1
                labelStr{results.bestFaciesOrder(j,k)} = data.faciesNames{k}; 
                colourCode(results.bestFaciesOrder(j,k),1:3) = data.faciesColours(k,2:4);
            end
            hX.XTickLabel = labelStr;
            hX.YTickLabel = labelStr;
            str = sprintf('Best M-value TP Matrix\nm=%4.3f', mStat);
            title(str);

            % Plot a best permutation facies coding order
            hY=subplot('Position',[0.85 0.10 0.1 0.8]);
            for k=1:maxNumbOfFacies
                yco = [k k k+1 k+1];
                xco = [1 2 2 1];
                faciesCol = [colourCode(k,1) colourCode(k,2) colourCode(k,3)];
                patch(xco, yco, faciesCol,'EdgeColor','black');
                text(1.25, k+0.2, labelStr(k), 'VerticalAlignment','bottom', 'HorizontalAlignment','left', 'FontSize',8);
                axis off;
            end
        end
    end
    
     % Complete analysis with a textbox message summary
    message = sprintf('For section %s\nCalculated %d facies coding premutations\nLargest markov stat value %5.4f found from %d permutations \nHighest m value * one-offset diagonal value %5.4f found from %d permutations\nFor optimised facies coding Markov order metric P Value %5.4f\n',...
        data.sectionName, results.permsCount, results.maxMarkov, results.maxMarkovCount, results.maxMarkovDiagProduct, results.maxDiagCount, markovPValueSum);
    m1 = msgbox(message,'Optimised cycle analysis');
end

% function [markovOrderMetricTest, oneOffsetDiagCheck, markovDiagProduct, maxMarkov, minMarkov, maxDiag, maxMarkovDiagProduct] = calcStatsForEachFaciesOrderPermutation(data, permsCount, sectionLength, allFaciesCombos)
function results = calcStatsForEachFaciesOrderPermutation(data, sectionLength, results)
    
    % Create a matrix containing the facies codes in all possible orders
    % allFaciesCombos should contain all the possible ways the n facies classes read in can be positioned in the TP matrix
    % In other words, a permutation 5 1 2 3 4 means that facies 1 is coded as facies 5 so that it positioned in row 5 of the TP matrix for this permutation
    results.allFaciesCombos = perms(data.faciesNumbers);
    results.permsCount = length(results.allFaciesCombos);

    % So this function needs to make a strat section in testFaciesSect using the facies coding in each row of the matrix and calculate the markov statistic for each testFaciesSect
    % Define the necessary variables to make a test section and record markov results from each test section. Pass this structure to the calculation functions
    results.testFaciesSect = zeros(1, sectionLength); % the facies section input but will contain facies coding from a single permutation
    results.markovOrderMetricTest = zeros(1,results.permsCount); results.permsCount
    results.oneOffsetDiagCheck = zeros(1,results.permsCount);
    results.markovDiagProduct = zeros(1,results.permsCount);
    results.maxMarkov = 0.0;
    results.minMarkov = 1000.0; % Max value is 1.0, so all values should be less than this, so good start value to find minimum
    results.maxDiag = 0.0;
    results.maxMarkovDiagProduct = 0.0;

    % Set variables to monitor the progress of the calculation
    progressCount = 0;
    progressIndicator = results.permsCount / 100.0;
    fprintf('Calculating M stats for %d permutations, * every %d iterations...\n', results.permsCount, progressIndicator);
    % Print a string of * symbols length permsCount as a marker for the * progress indicator printed below in the main calculation loop
    for j = 1:100.0
        fprintf('*');
    end
    fprintf('Done\n');
    
    % Loop through all of the facies coding permutations, code the vertical succession for each, and find the minimum, maximum m statistics, plus the max m-diagonal product
    for j = 1:results.permsCount

        % Make a vertical succession with each facies occurrence coded according to the facies coding in permutation j, not the original coding
        for k =1:sectionLength
            oneFacies = data.faciesCodes(k);
            results.testFaciesSect(k) = results.allFaciesCombos(j, oneFacies); % Swap the original facies code from the data section with a substitute code from permutation j
        end

        % Calculate the markov statistic and the one offset diagonal for this permutation
        [results.markovOrderMetricTest(j), results.oneOffsetDiagCheck(j)] = calculateTPMatrixAndOrderMetric(results.testFaciesSect, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01]));
        results.markovDiagProduct(j) = results.markovOrderMetricTest(j) * results.oneOffsetDiagCheck(j);

        % if the result is the highest markov stat yet, remember it
        if results.markovOrderMetricTest(j) >= results.maxMarkov
            results.maxMarkov = results.markovOrderMetricTest(j);
        end

        if results.oneOffsetDiagCheck(j) >= results.maxDiag
            results.maxDiag = results.oneOffsetDiagCheck(j);
        end

        if results.markovDiagProduct(j) >= results.maxMarkovDiagProduct
            results.maxMarkovDiagProduct = results.markovDiagProduct(j);
        end

        % if the result is the lowest markov stat yet, remember it 
        if results.markovOrderMetricTest(j) <= results.minMarkov
            results.minMarkov = results.markovOrderMetricTest(j);
        end
        
        
        % Output a progress indicator - useful for calculations with more permutations so that the user can see how far the calculation has got
        progressCount = progressCount + 1;
        if progressCount == progressIndicator
            fprintf('*');
            progressCount = 0;
        end
    end
end

function results = findOptimumAndWorstPermutations(data, sectionLength, maxNumbOfFacies, results)

    ROUND_ERROR = 0.0000001;
    
    results.maxMarkovPerm = zeros(1,results.permsCount); % Records the array index position of each test section that gives a maximum mvalue
    results.minMarkovPerm = zeros(1,results.permsCount); % Records the array index position of each test section that gives a minimum mvalue
    results.maxDiagPerm = zeros(1,results.permsCount); % Records the array index position of each test section that gives a maximum diagonal value
    results.bestFaciesOrder = zeros(1, maxNumbOfFacies); % An array of the facies code permutations that give the maximum markov stat value
    results.worstFaciesOrder = zeros(1, maxNumbOfFacies); % An array of facies code permutations that give the minimum markov stat value
    results.bestDiagFaciesOrder = zeros(1, maxNumbOfFacies); % An array of facies code permutations that give the maximum one offset diagonal test value
    results.maxMarkovCount = 0; % The number of permutations that give the maximum markov value
    results.minMarkovCount = 0;  % The number of permutations that give the minimum markov value
    results.maxDiagCount = 0;  % The number of permutations that give the maximum one offset diagonal test value
    results.bestFaciesSect = zeros(1, sectionLength); % An array of the succession permutations that give the maximum markov stat value
    results.worstFaciesSect = zeros(1, sectionLength); % An array of the succession permutations that give the maximum markov stat value
    results.bestDiagFaciesSect = zeros(1, sectionLength); % An array of the succession permutations that give the maximum markov stat value

    for j = 1:results.permsCount

        % Record information about the facies order permutations that give the highest M values
        % NB Convert to single precision for the >= test to remove possible rounding error effects that can otherwise introduced small
        % differences in markovOrderMetricTest that are not "real"
        if single(results.markovOrderMetricTest(j)) >= single(results.maxMarkov)
            results.maxMarkovCount = results.maxMarkovCount + 1;
            results.maxMarkovPerm(results.maxMarkovCount) = j;
            results.bestFaciesOrder(results.maxMarkovCount, :) = results.allFaciesCombos(j,:); % record the facies code permutation

            % Make a vertical succession with each facies occurrence coded according to the facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = data.faciesCodes(k);
                results.bestFaciesSect(results.maxMarkovCount, k) = results.allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end

        % Record information about the facies order permutations that give the lowest M values
        if single(results.markovOrderMetricTest(j)) <= single(results.minMarkov)
            results.minMarkovCount = results.minMarkovCount + 1;
            results.minMarkovPerm(results.minMarkovCount) = j;
            results.worstFaciesOrder(results.minMarkovCount, :) = results.allFaciesCombos(j,:); % record the facies code permutation

            % Make a vertical succession with each facies occurrence coded according to the  facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = data.faciesCodes(k);
                results.worstFaciesSect(results.minMarkovCount, k) = results.allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end
        
        % Record information about the facies order permutations with both a high m statistic and a high one-offset diagonal test value as indicated by a high markovDiagProduct
        if results.markovDiagProduct(j) >= results.maxMarkovDiagProduct - ROUND_ERROR
            results.maxDiagCount = results.maxDiagCount + 1;
            results.maxDiagPerm(results.maxDiagCount) = j;
            results.bestDiagFaciesOrder(results.maxDiagCount, :) = results.allFaciesCombos(j,:); % record the facies code permutation, gives same results as optimalCycleAnalysis
            
            % Make a vertical succession with each facies occurrence coded
            % according to the  facies coding from the jth permutation
            for k =1:sectionLength
                oneFacies = data.faciesCodes(k);
                results.bestDiagFaciesSect(results.maxDiagCount, k) = results.allFaciesCombos(j, oneFacies); % Make a complete vertical succession with facies coding from current permutation
            end
        end
    end
end

function gui = movingWindowAnalysisLoadAndPlot(gui, data, dataFileNameAndPath)

    fileIn = fopen(dataFileNameAndPath,'r');
    if (fileIn < 0)
        fprintf('\n\nPROBLEM: file %s could not be found to read\n', dataFileNameAndPath);
    else
        % Read and check the file header parameters
        headerInput = fscanf(fileIn,'%d',5);
        testSectLength = headerInput(1);
        minWindowSize = headerInput(2); 
        windowSizeInc = headerInput(3); 
        maxWindowSize = headerInput(4);
        windowBaseInc = headerInput(5);
    end
    
    if fileIn > 0 && testSectLength ~= data.sectLength
        fprintf('Problem reading file header - wrong maximum window size %d for the section length %d\n', testSectLength, data.sectLength);
    else
    
        totalWindowCount = 0;
        winDex = 1;
        while ~feof(fileIn)
            oneLineInput = fscanf(fileIn,'%d%d%f%f',4);
            
            if ~isempty(oneLineInput)
                
                windowSizeRecord(winDex) = oneLineInput(1); % Record the input values for one window size and position
                windowBaseIndex(winDex) = oneLineInput(2);
                markovOrderMetric(winDex) = oneLineInput(3); % NB elements 3 and 4 may be Nan if no p-value could be calculated for specific window
                markovPValueSum(winDex) = oneLineInput(4);
                         
                % because p values may be missing for some window sizes, check before trying to store all four values from the input line
                if ~isnan(oneLineInput(3)) % Nan is used as dummy value when a p value could not be calculated for any given window
                    
                    % All values present for this line, so use input values to calculate various things needed for the diagram plotting
                    windowTopIndex = windowBaseIndex(winDex) + (windowSizeRecord(winDex)-1); % -1 because windowBase starts at 1 not zero
                    windowBaseElev(winDex) = sum(data.faciesThick(1:windowBaseIndex(winDex)));
                    windowTopElev = sum(data.faciesThick(1:windowTopIndex));
                    windowMidElev(winDex) = (windowBaseElev(winDex) + windowTopElev) / 2.0;
                    meanThickUnitInWindow(winDex) = mean(data.faciesThick(windowBaseIndex(winDex):windowTopIndex)); % Used in the plotting routine to determine plot size of data point
                end

                winDex = winDex + 1;
                totalWindowCount = totalWindowCount + 1;
            end
        end

        fclose(fileIn);
 
        movingWindowAnalysisPlot(gui, data, totalWindowCount, minWindowSize, windowSizeInc, maxWindowSize, windowSizeRecord, windowMidElev, meanThickUnitInWindow, markovPValueSum);
    end
end

function gui = movingWindowAnalysisCalculateAndPlot(gui, data, dataFileNameAndPath)

    % Assume that the data filename has a .txt or .dat end, so remove this,
    % add a movingWindow data file ID label, then put .txt end back
    modifiedDataFileNameAndPath = strrep(dataFileNameAndPath,'.','MWindowAnalysis.');

    % Calculate the basic stats on the dataFacies array
    maxNumbOfFacies = length(data.faciesColours);
    results.testFaciesSect = 0; % results is passed to the functions below and returned full of results, so create now 

    minWindowSize = 10;
    maxWindowSize = data.sectLength;
    windowSizeInc = 2;
    windowBaseIncrement = 2;
    maxIterations = 200; % Number of monte carolo iterations to do to calculate the shuffled m value PDF
    
    windowSizesCount = (maxWindowSize - minWindowSize) / windowSizeInc;
    estimateTotalIterations = 0;
    for windowSizeLoop = minWindowSize: windowSizeInc: maxWindowSize
        estimateTotalIterations = estimateTotalIterations + ((data.sectLength - windowSizeLoop + 1) / windowBaseIncrement);
    end
    msgString = sprintf('%d iterations, will take estimated %4.3f hours assuming 200 Mc simulations per iteration\nTerminate now if too long and adjust window size parameters', ...
        estimateTotalIterations, ((estimateTotalIterations * 0.5)/60)/60); % Assuming average 0.5 seconds per iteration with 200 MC realizations
    msgbox(msgString);
    
    markovPValueSum = zeros(1,1);
    totalSectionThickness = sum(data.faciesThick(1:data.sectLength));
    
    fileOut = fopen(modifiedDataFileNameAndPath,'w');
    if (fileOut < 0)
        fprintf('\n\nWARNING: file %s could not be created, so no record will be made of this analysis\n', fName);
    end
    
    % Output the file 'header' information
    fprintf(fileOut,'%d\n%d %d %d\n%d\n',data.sectLength, minWindowSize, windowSizeInc, maxWindowSize, windowBaseIncrement);
    
    % Window order calculation should use the optimised facies order for this succession, therefore need to find this optimal order before proceeding
    % Calculate the various statistics returned from the function for each facies code permutation
    results = calcStatsForEachFaciesOrderPermutation(data, data.sectLength, results); 
   
    % now loop again to find and record all of the permutations and resulting strat sections that gave the maximum or minimum m and diag values found above
    % and record various information about each permutation
    results = findOptimumAndWorstPermutations(data, data.sectLength, maxNumbOfFacies, results);
    
    sectFaciesCodes = results.bestFaciesSect(1,:); % Set the facies code order used here to the first of the optimally-coded permutations
    
    totalWindowCount = 0;
    winDex = 1;
    
    for windowSize =  minWindowSize:windowSizeInc:maxWindowSize
        
        tic; % Time the calculations for each window size
        
        fprintf('Window size %d (max %d) starts unit 1 ends unit %d in increments of %d ', windowSize, maxWindowSize, (data.sectLength - windowSize) + 1, windowBaseIncrement);
        
        windowBaseIndex = 1; % All size windows have an initial base position at section unit 1
        oneWindowSizeCount = 0;
        while windowBaseIndex <= (data.sectLength - windowSize) + 1
            
            windowSizeRecord(winDex) = windowSize;
            windowTopIndex = windowBaseIndex + (windowSize-1); % -1 because windowBase starts at 1 not zero
            windowBaseElev(winDex) = sum(data.faciesThick(1:windowBaseIndex));
            windowTopElev(winDex) = sum(data.faciesThick(1:windowTopIndex));
            windowMidElev(winDex) = (windowBaseElev(winDex) + windowTopElev(winDex)) / 2.0;
            meanThickUnitInWindow(winDex) = mean(data.faciesThick(windowBaseIndex:windowTopIndex)); % Used in the plotting routine to determine plot size of data point
            
            windowSectFaciesCodes = sectFaciesCodes(windowBaseIndex: windowTopIndex); % Extract just the facies succession in the current window
            windowSectThicknesses = data.faciesThick(windowBaseIndex: windowTopIndex); % Not used except as dummy passes to the metric calculation functions
            numberOfSwaps = length(windowSectFaciesCodes); % e.g. if the window contains 10 facies units, 10 swaps ensures they are well-shuffled
            
            discreteFaciesCount = unique(windowSectFaciesCodes); % How many discrete facies in the window section? Need > 3 to calculate a meaningful m value
            
            if length(discreteFaciesCount) >= 4
             
                [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(windowSectFaciesCodes, data.faciesNames, 0); % 0 is plot flag, do not plot TP matrix results

                % Now shuffle the observed section and calculate the order metric and DW stat each time
                for j = 1:maxIterations
                    [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(windowSectFaciesCodes, windowSectThicknesses, numberOfSwaps);
                    [multiMarkovOrderMetric(j), oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(shuffledFacies, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01])); % note the zero is the TP matrix plot flag, set to false so no plot
                end

                % calculate P value from the shuffled section random models
                bins = 0:0.02:1.00; % because 0<m<=1. NB this gives 50 bins in the histogram array
                multiMarkovOrderMetricHistData = histc(multiMarkovOrderMetric, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq
        
                if (markovOrderMetric*50) >= 0.5 % Again 50 is the number of bins in the bar chart so multiply m by 50 to get array index value for frequencies
                    markovPValueSum(winDex) = sum(multiMarkovOrderMetricHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1 
                else
                    markovPValueSum(winDex) = 1.0; % Because value <0.5  will lead to array subscript error in line above, and value must be 1 since all PDF curve is at x>=0 on the histogram
                end

                if (fileOut > 0)
                    fprintf(fileOut,'%d %d %6.5f %6.5f\n', windowSize, windowBaseIndex, markovOrderMetric, markovPValueSum(winDex));
                end
            else
                if (fileOut > 0)
                    fprintf(fileOut,'%d %d NaN NaN\n', windowSize, windowBaseIndex);
                end
            end
            
            windowBaseIndex = windowBaseIndex + windowBaseIncrement;
            oneWindowSizeCount = oneWindowSizeCount + 1;
            totalWindowCount = totalWindowCount + 1;
            winDex = winDex + 1;
        end
        
        elapsedTime = toc;
        fprintf('done in %3.2f seconds (%5.4f seconds per window position)\n', elapsedTime, elapsedTime / double(oneWindowSizeCount));
    end
    
    fclose(fileOut);
   
    movingWindowAnalysisPlot(gui, data, totalWindowCount, minWindowSize, windowSizeInc, maxWindowSize, windowSizeRecord, windowMidElev, meanThickUnitInWindow, markovPValueSum);
end

function [meanPValuePerWindowSize, windowsEachSizeN] = calculateMeanPValuePerWindowSize(data, totalWindowCount, windowsEachSizeN, windowSizeRecord, markovPValueSum)
    
    % create a pValueMatrix sized based on range of window sizes, containing NaN values  as a starting condition
    pValueMatrix = nan(windowsEachSizeN, data.sectLength - windowSizeRecord(1));
    windowSize = windowSizeRecord(1); % Set windows size to the smallest window to start
    j = 1; % j is the array index for the window size, smallest window size j=1, largest window size j = windowsEachSizeN
    k = 1; % m is the vertical position of the window based in the section
    
    for m = 1:totalWindowCount
        pValueMatrix(j, k) = markovPValueSum(m);
        
        if (m+1 <= totalWindowCount && windowSizeRecord(m+1) ~= windowSize)
            windowSize = windowSizeRecord(m+1);
            j = j + 1; % increment j to the next window size
            k = 1; % Set k to position the window back at the base of the vertical succession
        else
            k = k + 1; % increment the vertical window position, NB +1 increment only, perhaps need to code actual window size increment?
        end
    end

    meanPValuePerWindowSize = zeros(1,windowsEachSizeN);
    for j=1:windowsEachSizeN
        oneWindowSizeTemp = pValueMatrix(j,:);
        oneWindowSizePValues = oneWindowSizeTemp(~isnan(oneWindowSizeTemp));
        meanPValuePerWindowSize(j) = mean(oneWindowSizePValues); % Calculate the mean p value for each window size, ignoring NaN blank values
    end
    
    fprintf('\n');
end

function movingWindowAnalysisPlot(gui, data, totalWindowCount, windowSizeMin, windowSizeInc, windowSizeMax, windowSizeRecord, windowMidElev, meanThickUnitInWindow, markovPValueSum)
% function movingWindowAnalysisPlot(gui, data, totalWindowCount, windowSizeMin, windowSizeInc, windowSizeRecord, windowMidElev, meanThickUnitInWindow, markovPValueSum)
    
    fprintf('Plotting moving window analysis results ');
    
    windowsEachSizeN = length(unique(windowSizeRecord)); 
    [meanPValuePerWindowSize, windowsEachSizeN] = calculateMeanPValuePerWindowSize(data, totalWindowCount, windowsEachSizeN, windowSizeRecord, markovPValueSum);

    % plot the vertical section on the left of the screen as a reference for the p values plot which will be on the right
    gui.movWindSectionPlot = subplot('Position', [0.04 0.1 0.075 0.85]);
    cla(gui.movWindSectionPlot);
    [gui,~] = plotObservedSection(gui, data, 1, data.sectLength);
    
    % Now plot the p values for each window size and position. Smallest windows on the left, largest on the right
    % Each window position is represented by one colour-coded p value square
    gui.movWindPValuesPlot = subplot('Position', [0.15 0.1 0.5 0.85]); % Plot each window as a point at the window midpoint elevation
    cla(gui.movWindPValuesPlot);
    windowPointXSize = windowSizeInc / 2;

    timerCount = 0; % because the plot may involve a large number of patches it might be slow, so show progress with a timer variable
    
    for x=1:totalWindowCount % Need to loop to plot a colour patch for each window size and position
        
        % Colour code the P values as follows:
        % p >= 0.10 no evidence for order - Red
        % 0.10 > p > 0.010 weak evidence for order - Yellow
        % 0.01 > p > 0.00 strong evidence for order - Green
        if markovPValueSum(x) < 0.0 || markovPValueSum(x) > 1.000001
            colour = [1 1 1];
        else if markovPValueSum(x) >= 0.10
                colour = [1.00 0.20 0.20];
            elseif markovPValueSum(x) >= 0.01
                colour = [1.00 1.00 0.20];
            else
                colour = [0.50 1.0 0.00];
            end
        end

        windowPointYSize = meanThickUnitInWindow(x) / 2; % Scale window size by the thickness of the strata within the window
        xco = [windowSizeRecord(x) - windowPointXSize, windowSizeRecord(x) - windowPointXSize, windowSizeRecord(x) + windowPointXSize, windowSizeRecord(x) + windowPointXSize];
        yco = [windowMidElev(x)-windowPointYSize, windowMidElev(x)+ windowPointYSize, windowMidElev(x)+ windowPointYSize, windowMidElev(x)- windowPointYSize];

        patch(xco, yco, colour, 'EdgeColor', 'none');
        
        timerCount = timerCount + 1; % Print a timer increment if we are another 10% through the whole plot
        if timerCount >= totalWindowCount / 10
            fprintf('*');
            timerCount = 0;
        end
    end
    
    % Want to also plot a mean p value for each window size, as a colour-coded horizontal bar along the top of the plot
    % Mean values are in meanPValuePerWindowSize so calculate a red-to-green colour from these values and plot it
    totatSectionThick = sum(data.faciesThick(:));
    meanValPlotBarThick = sum(data.faciesThick(int32(data.sectLength * 0.95):data.sectLength));
    
    for x = 1:windowsEachSizeN % So loop through each window size
        
        % Calculate the colours, including white if the mean value is not within the 0-1 required range
        if meanPValuePerWindowSize(x) < 0.0 || meanPValuePerWindowSize(x) > 1.000001
            colour = [1 1 1];
        else if meanPValuePerWindowSize(x) >= 0.10
                colour = [1.00 0.20 0.20];
            elseif meanPValuePerWindowSize(x) >= 0.01
                colour = [1.00 1.00 0.20];
            else
                colour = [0.50 1.0 0.00];
            end
        end
        
        % Coordinates for the horiztonal bar along the top of the main p value plot, each mean value in the correct position for its window size
        xco = [(windowSizeMin + ((x-1)*windowSizeInc)) - windowPointXSize, (windowSizeMin + ((x-1)*windowSizeInc)) - windowPointXSize, ...
            (windowSizeMin + (x*windowSizeInc)) - windowPointXSize, (windowSizeMin + (x*windowSizeInc)) - windowPointXSize] ;
        yco = [totatSectionThick, totatSectionThick+meanValPlotBarThick, totatSectionThick+meanValPlotBarThick, totatSectionThick];
        
        patch(xco, yco, colour, 'EdgeColor', 'none');
    end
    
    grid on;   
    xlabel('Window Size, number of total lithological units');
    
    %    Create a high-resultion 600dpi transparent background png file for this figure using export_fig which is an .m file in the Matlab folder  
    %   Usually commented out because it is slow and not necessary to do everytime code is run
%     export_fig MWAnalysisTrianglePlot.png -r600 -transparent
    
    % ScreenSize is a four-element vector: [left, bottom, width, height]:
    scrsz = get(0,'ScreenSize'); % vector 
    gui.movWindPValuesWindowSizeCrossPlot = figure('Position', [scrsz(3)*0.25, scrsz(4)*0.25, scrsz(3)*0.4, scrsz(4)*0.4]); 
    cla(gui.movWindPValuesWindowSizeCrossPlot);
    
    % Plot p value colour code zones
    xco = [0, 0, windowSizeMax, windowSizeMax]; % Could make the left margin windowSizeMin?
    ycoRed = [0.10, max(meanPValuePerWindowSize), max(meanPValuePerWindowSize), 0.10];
    ycoYellow = [0.01, 0.10, 0.10, 0.01];
    ycoGreen = [0.00, 0.01, 0.01, 0.00];
    patch(xco, ycoRed, [1.00 0.20 0.20], 'EdgeColor', 'none');
    patch(xco, ycoYellow, [1.00 1.00 0.20], 'EdgeColor', 'none');
    patch(xco, ycoGreen, [0.50 1.0 0.00], 'EdgeColor', 'none');
    
    for x = 1:windowsEachSizeN-1
        xco = [windowSizeMin + ((x - 1)*windowSizeInc), windowSizeMin + (x * windowSizeInc)];
        yco = [meanPValuePerWindowSize(x), meanPValuePerWindowSize(x+1)];
        
        line(xco, yco, 'color','black', 'LineWidth', 2);
%          line(xco, yco, 'color','black', 'LineWidth', 2, 'Marker', 'o'); % In most cases marker density is too high to be helpful
    end
    
    grid on;   
    ylim([-0.1 max(meanPValuePerWindowSize)+0.1])
    xlabel('Window Size (number of total lithological units in window)');
    ylabel('Mean p value');
    
    %    Create a high-resultion 600dpi transparent background png file for this figure using export_fig which is an .m file in the Matlab folder  
    %   Usually commented out because it is slow and not necessary to do everytime code is run
%     export_fig MWAnalysisPvsWindowSizePlot.png -r600 -transparent
    
    fprintf('all done\n');
end

function gui = calculateAndPlotStatsForOneWindow(gui, data, windowSize, windowPos)

    maxIterations = 200; % Number of monte carolo iterations to do to calculate the shuffled m value PDF
    maxNumbOfFacies = length(data.faciesColours);
    multiMarkovOrderMetric = zeros(1, maxIterations);
    results.testFaciesSect = 0; % results is passed to the functions below and returned full of results, so create now 
    
    % Window order calculation should use the optimised facies order for this succession, therefore need to find this optimal order before proceeding
    % Calculate the various statistics returned from the function for each facies code permutation
    results = calcStatsForEachFaciesOrderPermutation(data, data.sectLength, results); 
   
    % now loop again to find and record all of the permutations and resulting strat sections that gave the maximum or minimum m and diag values found above
    % and record various information about each permutation
    results = findOptimumAndWorstPermutations(data, data.sectLength, maxNumbOfFacies, results);
    
    sectFaciesCodes = results.bestFaciesSect(1,:); % Set the facies code order used here to the first of the optimally-coded permutations
    windowSectFaciesCodes = sectFaciesCodes(windowPos: (windowPos + windowSize)); % Extract just the facies succession in the current window

%     windowSectFaciesCodes = data.faciesCodes(windowPos: (windowPos + windowSize)); % Extract just the facies succession in the current window
    
    
    windowSectThicknesses = data.faciesThick(windowPos: (windowPos + windowSize)); % Not used except as dummy passes to the metric calculation functions
    numberOfSwaps = length(windowSectFaciesCodes); % e.g. if the window contains 10 facies units, 10 swaps ensures they are well-shuffled     
    discreteFaciesCount = unique(windowSectFaciesCodes); % How many discrete facies in the window section? Need > 3 to calculate a meaningful m value

    % Plot the vertical section from the specified window interval
    gui.sp1 = subplot('Position',[0.04 0.1 0.075 0.85]);
    cla(gui.sp1); % Clears the axis of whatever has been previously plotted
    [gui,~] = plotObservedSection(gui, data, windowPos, windowSize);
    
    % Calculate and plot the TP matrix for the strata in the window
    gui.sp2 = subplot('Position',[0.2 0.55 0.3 0.4]);
    cla(gui.sp2); % Clears the axis of whatever has been previously plotted
    [markovOrderMetric, oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(windowSectFaciesCodes, data.faciesNames, 1, gui.sp2); % 1 is plot flag, so plot the results
    
    % Calculate and plot the m statistic and the associated Monte Carlo test p value, but only if more than 3 distinct facies codes in the window section   
    if length(discreteFaciesCount) >= 4
           
        % Now shuffle the observed section and calculate the order metric and DW stat each time
        for j = 1:maxIterations
            [shuffledFacies, shuffledThick] = shuffleSectionNoSameTransitions(windowSectFaciesCodes, windowSectThicknesses, numberOfSwaps);
            [multiMarkovOrderMetric(j), oneOffsetMValueDiagCheck] = calculateTPMatrixAndOrderMetric(shuffledFacies, data.faciesNames, 0, subplot('Position',[0 0 0.01 0.01])); % note the zero is the TP matrix plot flag, set to false so no plot
        end

        % calculate P value from the shuffled section random models
        bins = 0:0.02:1.00; % because 0<m<=1. NB this gives 50 bins in the histogram array
        multiMarkovOrderMetricHistData = histc(multiMarkovOrderMetric, bins) / maxIterations; % Calculate frequency bins with histc but / by iterations to give relative freq

        if (markovOrderMetric*50) >= 0.5 % Again 50 is the number of bins in the bar chart so multiply m by 50 to get array index value for frequencies
            markovPValue = sum(multiMarkovOrderMetricHistData(int16(markovOrderMetric*50:length(bins)))); % area under curve from m to max m value 1 
        else
            markovPValue = 1.0; % Because value <0.5  will lead to array subscript error in line above, and value must be 1 since all PDF curve is at x>=0 on the histogram
        end
        
        plotMCMarkovHistogram(multiMarkovOrderMetricHistData, markovOrderMetric, [0.20 0.07 0.40 0.38]);
        
         % Complete analysis with a textbox message summary
        message = sprintf('Markov metric %4.3f\nMarkov order metric P Value %5.4f', markovOrderMetric,  markovPValue);
        m1 = msgbox(message,'Random model comparison');
    else
        msgbox('Too few unique facies in the window, cannot shuffle so cannot calculate p value\n');
    end
end

function gui = calculateSpectralAnalysis(gui, data)

    onePowerSpectrumAnalysis(data.faciesThick, length(data.faciesThick), 0, length(data.faciesThick));
end

function [glob] = onePowerSpectrumAnalysis(section, sectLength, plotFlag, windowSize)
% Calculates fft of thickness succession and test for significance of spectra peaks using an MC method
% We need to know though if the powers in this spectrum are significant.
% We can determine this using Monte Carlo analysis in which we run n iterations (e.g. n=1000), for each iterations shuffle the strata, then calculate a power spectrum
% for the shuffled strata. A data structure to store all these would be  big (e.g. 200*200*100*1000) so we do not want to store these. Instead calculate for each
% point, calculate a p value from the MC distribution and store this, since this is the final result we are interested in.
% We do want to be able to plot some examples though, so give the x-y coords of cases to plot, and when we come to them in the loop, plot them
    
    if sectLength > 0
        
%         significanceLimit = 0.01; % Peaks above this line have a very low (less than 1/100) probability of occurring by chance, and are therefore significant.
        iterationsMC = 1000;
        fprintf('Calculating %d Monte Carlo iterations...', iterationsMC);
        
        % Subtract mean thickness so that section mean is zero and data varies around mean - avoids nasty low-freq artefacts in power spectrum
        section = section - mean(section); 
        
        spectrumLength = round(sectLength / 2.0);

        % Set the various parameters and arrays required for the Monte Carlo analysis of the power spectra significance
        n = length(section);
        freqVects = NaN(1,n);
        powerSpectra = NaN(1,n);
        maxBins = 50; % The number of frequencies to be checked in power spectra
        powerSpectraMC = zeros(spectrumLength, iterationsMC); % size the power spectra arrays assuming a maximum power spectra length
        freqCountsMC = zeros(spectrumLength, maxBins);
        relativeFreqsMC  = zeros(spectrumLength, maxBins);
        binEdges = zeros(spectrumLength, maxBins);
        MCPDFBelowSpectrumPoint = zeros(1,maxBins);
        pValue = NaN(1,n);

        [sectPowerSpectrum, oneFreqVect] = pmtm(section, 4, sectLength, 1);

        % Get Monte Carlo iterations of the power spectrum for shuffled sections, to make a random model PDF
        powerSpectraMC(1:spectrumLength, 1:iterationsMC) = calculateShuffledMCPowerSpectra(section, sectLength, iterationsMC);

        freqCountsMC = zeros(spectrumLength, maxBins);

        % Loop through frequencies and count significant peaks at that frequency, and the p=0 point in the MC distribution of peaks
        for j = 1:spectrumLength
            % Need to put all MC realization values at each signal frequency j into a vector to calculate stats relative frequency distribution of those powers
            [tempFreqCounts, tempBinEdges] = histcounts(powerSpectraMC(j,1:iterationsMC), maxBins); % define constant number of bins for all frequencies - can then plot easily
            freqCountsMC(j,1:length(tempFreqCounts)) = tempFreqCounts;
            binEdges(j,1:length(tempBinEdges)) = tempBinEdges;
            relativeFreqsMC(j, :) = freqCountsMC(j, :) / sum(freqCountsMC(j,:)); % Convert frequency to relative frequency

            MCPDFBelowSpectrumPoint = (binEdges(j,1:length(tempBinEdges)-1) < sectPowerSpectrum(j)); % Should return a vector of 1s for bin edge values < jth spectrum peak
            pValue(j) = 1 - sum(MCPDFBelowSpectrumPoint .* relativeFreqsMC(j,:));
        end
    else     
        fprintf("Cannot analyse a zero-length section\n");
    end
    
    plotPowerSpectrumAndMCAnalysis(oneFreqVect, sectPowerSpectrum, freqCountsMC, binEdges, relativeFreqsMC, pValue);
end

function  powerSpectraMC = calculateShuffledMCPowerSpectra(section, standardSectLength, iterationsMC)
% Calculate power spectra for iterationsMC shuffled versions for the vertical section passed to this function
% This is then a Monte Carlo model of the power spectra, useful as a random model to indicate which peaks are significant in the spectrum

    spectrumLength = round(standardSectLength / 2.0);
    powerSpectraMC = zeros(1, iterationsMC);
    for m = 1:iterationsMC
        fprintf('%4d',m);
        section = section(section~=0); % Remove all the zero values, including padding zeros used to increase the section length to the standard length
        nThicknesses = length(section); %number of datapoints in series
        if nThicknesses > 0 % ensure shortest section are ecluded - check threshold length is same as calling functions
            shuffledSection = shuffleSection(section, standardSectLength); % note section passed with zero-thickness values removed, variable legnth and no padding
            [onePowerSpectrum, ~] = pmtm(shuffledSection, 4, length(shuffledSection), 1);
            powerSpectraMC(1:spectrumLength, m) = onePowerSpectrum(1:spectrumLength);
        else
            powerSpectraMC(1:spectrumLength, m) = 0;
        end
        fprintf('\b\b\b\b');
    end
    fprintf('Done\n');
end

function shuffledSection = shuffleSection(section, standardSectLength)
% Shuffle the facies succession to ensure a random configuration, and pad
% with zeros after shuflling (because we don't want to shuffle a padded
% section because it would introduce lots of zeros into the section

    % Shuffle the section using randperm - ~ one order of magnitude faster than individual element swaps in a for loop!
    shuffledSection = section(randperm(numel(section)));
    shuffledSection = reshape(shuffledSection, [1,length(section)]);
    
    zeroPadding = zeros(1,(standardSectLength - length(shuffledSection)));
    shuffledSection = horzcat(shuffledSection, zeroPadding);
end

function plotPowerSpectrumAndMCAnalysis(oneFreqVect, onePowerSpectrum, freqCountsMC, binEdges, relativeFreqsMC, pValues)

    scrsz = get(0,'ScreenSize'); % screen dimensions vector
    ffPSpect = figure('Visible','on','Position',[100, 0, scrsz(3)*0.5, scrsz(4)*0.95]);
    
    subplot(2,1,1);
    hold on
    
    oneLayersVect = 1 ./ oneFreqVect; % Convert frequency to number of layers for plotting
    plotLimitNLayers = 100; % Want to analyse cyclicity at 100 layers wavelength or less
    startPos = numel(find(oneLayersVect >= plotLimitNLayers )); % Find the first position in the vector where layers<=100
    oneLayersVectPlot = oneLayersVect(startPos:numel(oneLayersVect));
    onePowerSpectrumPlot = onePowerSpectrum(startPos: numel(oneLayersVect));
    pValues = pValues(startPos: numel(oneLayersVect));
    
    maxRelFreqMC = 0.1; 
%     numberOfPValueLabels = 3;
%     maxPowerOnPlot = max(max(binEdges));
 

    % Now plot the relative frequencies of powers for each frequency produced by the MC analysis of the random shuffled strata
    for j = 2: length(oneLayersVectPlot)-1 % Loop through all frequencies included in the power spectrum, so the x axis of the power spectrum plot
        freqDimsMC = size(freqCountsMC(:,:));
        halfXBinSize = (oneLayersVectPlot(j-1) - oneLayersVectPlot(j+1)) / 2.0;
        xcoPlot = [oneLayersVectPlot(j) - halfXBinSize, oneLayersVectPlot(j) - halfXBinSize, oneLayersVectPlot(j) + halfXBinSize, oneLayersVectPlot(j) + halfXBinSize];
            
        for k = 1: freqDimsMC(2) % Loop through the bins created for the PDF, from the size of the 2nd dimension of freqDimsMC
            ycoPlot = [binEdges(j,k), binEdges(j,k+1), binEdges(j,k+1), binEdges(j,k)];

            if freqCountsMC(j,k) > 0
                if relativeFreqsMC(j,k) > maxRelFreqMC
                    redPinkScale = 0;
                else
                    redPinkScale = 1-(relativeFreqsMC(j,k)/maxRelFreqMC);
                end
                patch(xcoPlot, ycoPlot, [1 redPinkScale, redPinkScale*0.6] , 'EdgeColor','none');
             end
        end
    end
    
    % Plot the line showing the powers across the range of layer numbers
    line(oneLayersVectPlot, onePowerSpectrumPlot, 'Color', [0, 0, 0], 'LineWidth',1);
    
%     % circle and label P values for the most significant power spectra peaks
%     [powerVals, powerValIndices] = maxk(onePowerSpectrumPlot, numberOfPValueLabels);
%     scatter(oneLayersVectPlot(powerValIndices), powerVals, 'o','b')
%     for j = 1:numberOfPValueLabels
%         txt = sprintf("p=%4.3f", pValues(powerValIndices(j)));
%         text(oneLayersVectPlot(powerValIndices(j)), powerVals(j) + maxPowerOnPlot / 50, txt)
%     end

    ax = gca;
    ax.FontSize = 14;
    xlabel('Number of strat layers');
    ylabel('Spectral power');
    grid on;
    hold off
    
    % Now plot the p-values for the same range of number of layers
    subplot(2,1,2);
    
    % plot a green vertical rectangle over p values < 0.01
    for j = 1:numel(oneLayersVectPlot)-1
        
        if pValues(j) < 0.01 && pValues(j+1) < 0.01
            colourVector = [0, 1, 0];
            fprintf("p=%5.4f at beds count frequency %3.2f is strong evidence for cyclicity on this scale\n", pValues(j), oneLayersVectPlot(j));
        elseif pValues(j) < 0.1 && pValues(j+1) < 0.1
            colourVector = [1, 1, 0];
        else
            colourVector = [1, 0, 0];
        end
       
       xco = [oneLayersVectPlot(j), oneLayersVectPlot(j+1), oneLayersVectPlot(j+1), oneLayersVectPlot(j)];
       yco = [0,0,1,1];
       patch(xco, yco, colourVector, "FaceAlpha", 0.25, "EdgeColor", "none");
    end
    hold on
    
    % plot the p values as a thick black line
    line(oneLayersVectPlot, pValues, 'Color', [0, 0, 0], 'LineWidth',1);
    
    xlabel('Number of strat layers');
    ylabel('p value');
    grid on;
    ax = gca;
    ax.FontSize = 14;
end
