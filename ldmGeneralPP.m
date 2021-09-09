function ldmGeneralPP(desc,trncTag,trncLim,cutTag,cutLim,slcTag,slcVal)
% LDMGENERALPP(desc,trncTag,trncLim,cutTag,cutLim,slcTag,slcVal)
%
% Post-processes the fully 3D datasets (density, melt fraction, and 
% temperature) generated from the Laguna del Maule magmatic simulations and
% forward modeling (electrical resistivity and shear velocity) thereof.
% Post-processing of forward modeled gravity anomaly is done elsewhere.
%
% For every timestep, a PNG figure representing the full 3D dataset as a 
% series of preset isosurfaces is produced. Per the input provided, 
% accessory figures such as vertical truncations, interior cutaways, and 
% 2D-slices of the data may also be produced. All timesteps of each figure
% set are compiled into AVI-format movies. All files are labeled with 
% two-part ID tags denoting (1) the subset of data depicted (e.g., the 
% "Full" dataset or a 2D-slice); and (2) Any cosmetic  modifications (e.g.,
% truncation or cutaways). These partial tags are provided via function
% input along with the accessory figure specifications.
%
% NOTE: The 3D axes convention of the original magmatic simulations is
% maintained, meaning that the Y-axis is vertical (and negative from the
% surface) and the X- and Z-axes are the two lateral axes.
%
% INPUT:
% 
% desc      Indicates the data to be post-processed. Shorthand options are: 
%              'drho'   for Density [kg/m3]
%              'erho'   for Electrical Resistivity [Ohm m]
%              'f'      for Melt Fraction [-]
%              'vp'     for P-wave velocity [km/s]
%              'vs'     for S-wave Velocity [km/s]
%              't'      for Temperature [degree C]
%              'deltat' for Temperature Anomaly [degree C]
% trncTag   ID tag labeling vertical truncation accessory figures
%              Default: {'Trnc'}
% trncLim   Y-limits [km] for each requested truncation 
%              Default: [-20 0] 
% cutTag    ID tag labeling interior cutaway accessory figures
%              Default: {'CutX' 'CutZ'}
% cutLim    X,Y,Z-limits [km] for each requested cutaway 
%              Default: [15 35 -20 0 25 35;
%                        25 35 -20 0 15 35]
% slcTag    ID tag labeling 2D slice accessory figures
%              Default: {'SlcX' 'SlcZ' 'SlcY'}
% slcVal    X/Y/Z-value [km] for each requested slice. Only slices normal
%              to a principle axis are allowed (NaNs denote the other axes)
%              Default: [25 NaN NaN;
%                        NaN NaN 25;
%                        NaN -5 NaN]
%
% SEE ALSO: LDMINPUTPROCESSING
%
% Last modified by gleggers-at-gatech.edu, 10/18/2020

% INPUT VARIABLE INTERPRETATION
disp('Now beginning post-processing routine...')
tic

% Default values for input variables
% Truncation Figures: Vertical axis limits and assocated tags for labeling
defval('trncTag',{'Trnc'});
defval('trncLim',[-20 0]);            % Truncation to near-surface [km]
% Cutaway Figures: Axes limits and associated tags for labeling
defval('cutTag',{'CutX' 'CutZ'});
defval('cutLim',[15 35 -20 0 25 35;   % Cutaway parallel to X-axis [km]
                 25 35 -20 0 15 35]); % Cutaway parallel to Z-axis [km]
% 2D Slice Figures: Axis value for slice and associated tags for labeling
defval('slcTag',{'SlcX' 'SlcZ' 'SlcY'});
defval('slcVal',[25 NaN NaN;          % Slice along X = 25 km
                 NaN NaN 25;          % Slice along Z = 25 km
                 NaN -5 NaN]);        % Slice along Y = -5 km

% Note: Should do a quality check of inputs, but not for now.
             
% Other relevant variables
figv = 'off';        % Figure visibility switch
rndr = 'painters';  % Renderer
fres = '-r800';     % Figure writeout resolution

% READ IN AND PROCESS DATA
disp('...Processing inputs and reading in data')

% "Half" tag that indicates the full dataset was imaged
fllTag = 'Full';

% Conduct the standard LdM input processing
[dataAll,bnd,res,tStep,iVal,iSort,cref,tagDat,daName,daUnit,flux,reali] =...
    ldmInputProcessing(desc);

% Parcel out the physical boundaries of the simulation space
xbnd = bnd(1,:); ybnd = bnd(2,:); zbnd = bnd(3,:);
% Parcel out the data resolution
dx = res(1); dy = res(2); dz = res(3);

% Define grid vectors giving centerpoints of cells in the 3D domain
Xgv = (xbnd(1)+dx/2):dx:(xbnd(2)-dx/2);
Ygv = (ybnd(1)+dy/2):dy:(ybnd(2)-dy/2);
Zgv = (zbnd(1)+dz/2):dz:(zbnd(2)-dz/2);

% Define a meshed grid for that domain [Note: Y is the vertical dimension]
[Xg,Zg,Yg] = meshgrid(Xgv,Zgv,Ygv);

% Switch on a flag if processing magnetotelluric data
if contains(tagDat,'eResistivity')
    flagMT = 1;
else
    flagMT = 0;
end

% Prepare a lower-cased data name
daNameL = lower(daName);
% But recapitalize the wave letter for seismic velocities
if contains(tagDat,{'pVelocity' 'sVelocity'})
    daNameL = [upper(daNameL(1)) daNameL(2:end)];
end

% GENERATE AND SAVE POST-PROCESSING FIGURES
fprintf('...Conducting post-processing of %s data\n',daName)
disp('...Generating figures')

% Sort isovalues from what would be outermost to innermost surface
iVal = sort(iVal,iSort);

% Alpha gradient for isosurfaces from most transparent to opaque
iTran = linspace(0.2,1,length(iVal));

% Color axis range and colorbar ticks/labels for 2D accessory figures
cax = [min(iVal) max(iVal)];
cbTick = sort(iVal,'ascend'); cbTickL = num2str(cbTick(:));

% Tick mark spacing
dtickC = 10; % Coarser tick mark spacing [km] for the full section 
dtickF = 5;  % Finer tick mark spacing [km] for the accessory figures 

% Flag to signal first time through loop
flagLoop = 0;

% Visualize the shear velocity at each timestep
for i = 1:length(tStep)
    
    % Make a timestep ID tag
    tagTim = sprintf('T%03d',i);
    
    % Extract data for the timestep
    data = dataAll(:,:,:,i);
    
    % Start a figure counter for this timestep
    indFig = 1;
    
    % Open a blank figure for the 3D isosurfaces of the full dataset
    ffh = figure('Renderer',rndr,'Visible',figv);
    afh = gca;
    
    % For each isovalue...
    for j = 1:length(iVal)
        
        % Compute the isosurface geometry (faces/vertices) at that isovalue
        FV = isosurface(Xg,Zg,Yg,data(:,:,:),iVal(j));
        % Create filled polygons using those calculated faces and vertices
        p = patch(FV);
        % Calculate and set the normals to the isosurface vertices
        isonormals(Xg,Zg,Yg,data(:,:,:),p)
        % Set a color index and transparency for the isosurface
        set(p,'CData',j,'FaceAlpha',iTran(j));
        % Tweak isosurface for (1) uniform colors; (2) no edges displayed;
        % and (3) direct mapping between color indices and colormap
        set(p,'FaceColor','flat','EdgeColor','none',...
              'CDataMapping','direct');
    end

    % Get index of any closing parenthesis in the colormap reference string
    pind = strfind(cref,')');
    
    % Insert length of isovalues array into the string reference
    % For nonempty index (due to FLIP), inset before closing parenthesis
    if ~isempty(pind)
        iCref = sprintf('%s(%d)%s',cref(1:pind-1),length(iVal),...
                cref(pind:end));
    
    % Otherwise, take the length onto the end of the string reference
    else
        iCref = sprintf('%s(%d)',cref,length(iVal));
    end
    
    % Make a limited colormap with warm colors indicating low velocity
    cmap = eval(iCref);
    % For isosurfaces drawn in descending order...
    if strcmp(iSort,'descend')
        % Flip the colormap so it maps from high to low values
        cmap = flip(cmap); 
    end
    % Set the limited colormap with direct mapping to the color indices
    colormap(cmap);
    
    % Set order for rendering surfaces to the order they were drawn
    set(afh,'SortMethod','childorder')
    
    % Plot accoutrements
    % Set the aspect ratio to an equal [1 1 1]
    axis(afh,'equal')
    % Set axes limits to extent of modeled data (minding the 3D convention)
    xlim(afh,xbnd); ylim(afh,zbnd); zlim(afh,ybnd);
    % Set axes tick marks
    xticks(afh,xbnd(1):dtickC:xbnd(2));
    yticks(afh,zbnd(1):dtickC:zbnd(2));
    zticks(afh,ybnd(1):dtickC:ybnd(2));
    % Give axes labels
    xlabel(afh,'distance X (km)','FontSize',12);
    ylabel(afh,'distance Z (km)','FontSize',12);
    zlabel(afh,'depth Y (km)','FontSize',12);
    % Set the box and grid parameters
    box on; grid on;
    
    % Add a colorbar with labeled isovalues (and adjust up a half-step)
    cbh = colorbar(afh,'YTick',(1:length(iVal))+0.5,...
                       'YTickLabel',num2str(iVal(:)));
    % Add a colorbar label with our without units as appropriate
    if isempty(daUnit)
        ylabel(cbh,sprintf('%s',daNameL),'FontSize',12);
    else
        ylabel(cbh,sprintf('%s [%s]',daNameL,daUnit),'FontSize',12);
    end
    % For isosurfaces drawn in descending order...
    if strcmp(iSort,'descend')
        % Reverse the direction of the colorbar
        set(cbh,'Direction','reverse')
    end
    
    % Give a title
    tstr = sprintf('LdM (F=%02d, R=%03d) %s\nt=%06.1f ky',...
        flux,reali,daName,tStep(i));
    th = title(afh,tstr,'FontSize',14);
    
    % Tweak the viewing geometry and lighting
    view([-50,16]);             % Azimuth [deg] and elevation [deg]
    camproj orthographic        % Maintain actual object sizes and angles
    camlight; lighting gouraud; % Suitable for curved surfaces
    
    % Craft and retain a figure ID tag on the first timestep
    if flagLoop == 0
        tagFig{indFig} = fllTag;
    end

    % Name the full data figure with no modifications
    fname = sprintf('%s_%s_%s',tagDat,tagFig{indFig},tagTim);
    % Print the RGB image data with high resoultion
    cdata = print(ffh,'-RGBImage',fres);
    % Save the RGB data to file and retain it as a movie frame
    imwrite(cdata,sprintf('%s.png',fname),'PNG');
    F(indFig,i) = im2frame(cdata);
    
    % Increase the figure index
    indFig = indFig + 1;
    
    % Set a finer vertical tick mark spacing for the truncation figures
    zticks(afh,ybnd(1):dtickF:ybnd(2));
   
    % For each requested truncation accessory figure...
    for j = 1:length(trncTag)
        
        % Apply the truncation the full data figure
        zlim(afh,trncLim(j,:));
        
        % Craft and retain a figure ID tag on the first timestep
        if flagLoop == 0
            tagFig{indFig} = sprintf('%s%s',fllTag,trncTag{j});
        end
        
        % Name the modified full data figure
        fname = sprintf('%s_%s_%s',tagDat,tagFig{indFig},tagTim);
        % Print the RGB image data with high resoultion
        cdata = print(ffh,'-RGBImage',fres);
        % Save the RGB data to file and retain it as a movie frame
        imwrite(cdata,sprintf('%s.png',fname),'PNG');
        F(indFig,i) = im2frame(cdata);        
        
        % Increase the figure index
        indFig = indFig + 1;
    end
    
    % Set a finer lateral tick mark spacing for the cutaway figures
    xticks(afh,xbnd(1):dtickF:xbnd(2));
    yticks(afh,zbnd(1):dtickF:zbnd(2));
    
    % For each requested cutaway accessory figure...
    for j = 1:length(cutTag)
        
        % Apply the cutaway limits to the full data figure
        xlim(cutLim(j,1:2)); ylim(cutLim(j,5:6)); zlim(cutLim(j,3:4));
        
        % Craft and retain a figure ID tag on the first timestep
        if flagLoop == 0
            tagFig{indFig} = sprintf('%s%s',fllTag,cutTag{j});
        end
        
        % Name the modified full data figure
        fname = sprintf('%s_%s_%s',tagDat,tagFig{indFig},tagTim);
        % Print the RGB image data with high resoultion
        cdata = print(ffh,'-RGBImage',fres);
        % Save the RGB data to file and retain it as a movie frame
        imwrite(cdata,sprintf('%s.png',fname),'PNG');
        F(indFig,i) = im2frame(cdata);

        % Increase the figure index
        indFig = indFig + 1;
    end

    % Close the figure for the 3D isosurfaces
    close(ffh);
    
    % For each requested 2D slice accessory figure...
    for j = 1:length(slcTag)
      
        % Index of axis through which the perpendicular slice cuts
        indAx = find(~isnan(slcVal(j,:)));
        % Get the axis value through which the slice cuts
        slc = slcVal(j,indAx);

        % Based on which axis that is...
        switch indAx
            
            % Slice through the X-axis yielding a ZY-plane
            case 1
                
                % Reassign the slice values, changing NaNs as empty arrays
                xSlc = slcVal(j,indAx); ySlc = []; zSlc = [];
                % Order to which to permute the data slice dimensions
                pOrder = [2 1];
                
                % Horizontal/Vertical (for the 2D plot) data bounds
                hbnd = zbnd; vbnd = ybnd;
                % Horizontal/Vetical (for the 2D plot) data resolutions
                dh = dz; dv = dy;                
                % Title and horizontal/vertical axes label components
                titl = 'X'; haxl = 'distance Z'; vaxl = 'depth Y';
                
            % Slice through the Y-axis yielding an XZ-plane
            case 2
                
                % Reassign the slice values, changing NaNs as empty arrays
                xSlc = []; ySlc = slcVal(j,indAx); zSlc = [];
                % Order to which to permute the data slice dimensions
                pOrder = [1 2];
                
                % Horizontal/Vertical (for the 2D plot) data bounds
                hbnd = xbnd; vbnd = zbnd;
                % Horizontal/Vetical (for the 2D plot) data resolutions
                dh = dx; dv = dz;                
                % Title and horizontal/vertical axes label components
                titl = 'Y'; haxl = 'distance X'; vaxl = 'distance Z';

            % Slice through the Z-axis yielding an XY-plane
            case 3
                
                % Reassign the slice values, changing NaNs as empty arrays
                xSlc = []; ySlc = []; zSlc = slcVal(j,indAx);
                % Order to which to permute the data slice dimensions
                pOrder = [2 1];
                
                % Horizontal/Vertical (for the 2D plot) data bounds
                hbnd = xbnd; vbnd = ybnd;
                % Horizontal/Vetical (for the 2D plot) data resolutions
                dh = dx; dv = dy;                
                % Title and horizontal/vertical axes label components
                titl = 'Z'; haxl = 'distance X'; vaxl = 'depth Y';
        end

        % Open a dummy figure with its visibility off
        dfh = figure('Visible','off');
        
        % Generate the requested 2D slice of the 3D data
        sh = slice(Xg,Zg,Yg,data,xSlc,zSlc,ySlc);
        % Extract the slice color data, which is the actual 2D data
        dataSlc = sh.CData;         
        
        % Close the dummy figure
        close(dfh);
        
        % Permute the data slice dimensions
        dataSlc = permute(dataSlc,pOrder);
        
        % If processing magnetotelluric data...
        if flagMT == 1
            
            % Take a base 10 logarithm of the data for scaling
            dataSlc = log10(dataSlc);
            
            % If on the first timestep and first 2D slice...
            if (i == 1 && j == 1)
                
                % Take a log of color axis parameters
                cax = log10(cax);
                cbTick = log10(cbTick);                
            end
        end
        
        % Open a blank figure for the 2D slices of the full dataset
        fsh = figure('Visible',figv);
        ash = gca;
        
        % Image sliced shear velocity (with correct cell center locations)
        ph = imagesc(hbnd+[dh -dh]/2,vbnd+[dv -dv]/2,dataSlc,cax);
        
        % Set the colormap
        cmap = eval(cref);
        colormap(cmap);
        
        % Plot accoutrements
        % Set the vertical axis to have low values at the bottom
        set(ash,'Ydir','normal')
        % Set the aspect ratio to an equal [1 1 1]
        axis(ash,'equal')
        % Set axes limits to extent of sliced data
        xlim(ash,hbnd); ylim(ash,vbnd);
        % Set axes tick marks
        xticks(ash,hbnd(1):dtickC:hbnd(2));
        yticks(ash,vbnd(1):dtickC:vbnd(2));
        % Give axes labels
        xlabel(ash,sprintf('%s (km)',haxl),'FontSize',12);
        ylabel(ash,sprintf('%s (km)',vaxl),'FontSize',12);
        % Add a labeled colorbar
        cbh = colorbar(ash,'YTick',cbTick,'YTickLabel',cbTickL);
        % Add a colorbar label with our without units as appropriate
        if isempty(daUnit)
            ylabel(cbh,sprintf('%s',daNameL),'FontSize',12);
        else
            ylabel(cbh,sprintf('%s [%s]',daNameL,daUnit),'FontSize',12);
        end       
        % Give a title
        tstr = sprintf('LdM (F=%02d, R=%03d) %s\n%s=%02d km, t=%06.1f ky',...
            flux,reali,daName,titl,slc,tStep(i));
        th = title(ash,tstr,'FontSize',14);
        
        % Craft and retain a figure ID tag on the first timestep
        if flagLoop == 0
            tagFig{indFig} = slcTag{j};
        end
        
        % Name the sliced data accessory figure
        fname = sprintf('%s_%s_%s',tagDat,tagFig{indFig},tagTim);
        % Print the RGB image data with high resoultion
        cdata = print(fsh,'-RGBImage',fres);
        % Save the RGB data to file and retain it as a movie frame
        imwrite(cdata,sprintf('%s.png',fname),'PNG');
        F(indFig,i) = im2frame(cdata);

        % Increase the figure index
        indFig = indFig + 1;
        
        % Make truncated figures if the slice is through the X or Z-axes
        if (indAx == 1) || (indAx == 3)
        
            % Set a finer vertical tick mark spacing
            yticks(ash,ybnd(1):dtickF:ybnd(2));
            
            % For each requested truncation accessory figure...
            for k = 1:length(trncTag)
        
                % Apply the truncation the sliced data figure
                ylim(ash,trncLim(k,:));
        
                % Craft and retain a figure ID tag on the first timestep
                if flagLoop == 0 
                    tagFig{indFig} = sprintf('%s%s',slcTag{j},trncTag{k}); 
                end
                
                % Name the modified sliced data figure
                fname = sprintf('%s_%s_%s',tagDat,tagFig{indFig},tagTim);
                % Print the RGB image data with high resoultion
                cdata = print(fsh,'-RGBImage',fres);
                % Save the RGB data to file and retain it as a movie frame
                imwrite(cdata,sprintf('%s.png',fname),'PNG');
                F(indFig,i) = im2frame(cdata);

                % Increase the figure index (for later potential figures)
                indFig = indFig + 1;
            end
        end
        
        % Close the figure for the 2D slice
        close(fsh);
    end
    
    % Update the loop flag at the end of the first loop
    if flagLoop == 0
        flagLoop = 1;
    end
    
    % Give an update
    fprintf('......Finished timestep %d of %d\n',i,length(tStep))
end


% AGGREGATE FIGURES INTO MOVIES
disp('...Aggregating figures into movies')

% For every figure set generated...
for i = 1:length(tagFig)
    
    % Name the movie
    vname = sprintf('%s_%s',tagDat,tagFig{i});
    
    % Set the movie object properites
    mov = VideoWriter(vname,'Motion JPEG AVI');
    mov.FrameRate = 1; % Frames per second
    
    % Save the movie
    open(mov);
    writeVideo(mov,F(i,:));
    close(mov)
end


% END OF FUNCTION CLOSE OUT
tElapse = toc;
tElapseHrs = tElapse/3600;
fprintf('%s post-processing completed\n',daName)
fprintf('Elapsed time was %6.0f seconds, or %0.4f hours\n\n',...
    tElapse,tElapseHrs)
end

% AUXILIARY SUBFUNCTIONS

function defval(name,value)
% DEFVAL(name,value)
%
% Assigns a default value to the named variable
%
% INPUT:
% 
% name    A string, enclosed in single quotes, with a variable name
% value   The value, whatever it is, that you want the variable to have 
%
% OUTPUT:
%
%      None. The variables appear as if by magic into your workspace or
%      will be available inside your function.
%
% NOTE: 
%
% This won't work for an unassigned structure variable.
%
% Last modified by ebrevdo-at-alumni-princeton.edu, 05/28/2011
% Last modified by fjsimons-at-alum.mit.edu, 12/20/2012
% 
% It appears that defval('bla',functioncall) evaluates the function call
% regardless of whether or not 'bla' has been assigned.

if ~ischar(name),
  error(sprintf(['The first argument of DEFVAL ',...
		'has to be a string with a variable name']));
end

% Always do it is our default here
si=1;
% If it exists... as a variable (say it, it makes it faster!)
if evalin('caller',[ 'exist(''' name ''',''var'')']);
  % ... and it's empty, do it; but don't do it if it's non empty
  si=evalin('caller',[ 'isempty(' name ')']);
end
% Do it or not
if si
  assignin('caller',name,value);
end
end
