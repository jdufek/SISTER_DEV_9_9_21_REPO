function zrSatDiffPP
%
% Make figures of the zircon oversaturation for a Laguna del Maule
% simulation

% Read the INPUT file for the simulation
fid = fopen('INPUT');
% Simulation domain [cells]
xn = cell2mat(textscan(fid,'%d','Headerlines',1));
yn = cell2mat(textscan(fid,'%d','Headerlines',1));
zn = cell2mat(textscan(fid,'%d','Headerlines',1));
% Physical domain [m]
xdim = cell2mat(textscan(fid,'%d','Headerlines',1));
ydim = cell2mat(textscan(fid,'%d','Headerlines',1));
zdim = cell2mat(textscan(fid,'%d','Headerlines',1));
% Time parameters [yr]
dt = cell2mat(textscan(fid,'%d','Headerlines',3));
timesteps = cell2mat(textscan(fid,'%d','Headerlines',1));
twriteout = cell2mat(textscan(fid,'%d','Headerlines',1));
fclose(fid);

% Read in the start time
fid = fopen('start_time');
tstart = double(cell2mat(textscan(fid,'%d')));
fclose(fid);

% Number of written out timesteps
tn = timesteps/twriteout;

% Extract information from the directory name
[~,dName] = fileparts(pwd);
ind = strfind(dName,'_');
flux = str2double(dName(ind(2)+2:ind(3)-1)); % Flux of the simulation
reali = str2double(dName(ind(3)+2:end)); % Realization of the simulation

% Convert spatial domain
xdim = xdim/1e3; % [m --> km]
ydim = ydim/1e3; % [m --> km]
zdim = zdim/1e3; % [m --> km]

% Define grid vectors for a meshgrid (remember 3D orientation)
Xgv = [0.5:1:double(xn)-0.5]*double(xdim)/double(xn);
Ygv = [0.5:1:double(yn)-0.5]*double(ydim)/double(yn);
Zgv = [0.5:1:double(zn)-0.5]*double(zdim)/double(zn);
[X,Y,Z] = meshgrid(Xgv,Zgv,Ygv);

% Read in melt fraction data
% If a MATLAB archived version exists, read from that
if exist('zrSaturation.mat','file')
    load('zrSaturation.mat');

% Otherwise, read it in anew
else
    zrAll = load('ZrsatDif');
    
    % Save a MATLAB archived version for quicker retrieval
    save('zrSaturation.mat','zrAll')
end

% Reshape and permute it into the correct form for plotting
zrAll = reshape(zrAll,[xn yn zn tn]);
zrAll = permute(zrAll,[3 1 2 4]);

% Melt fraction isovalue surfaces to image
iSat = [0 10 20 30 40];
% Transparency gradient (0.2 to 1) for isosurfaces
iTran = linspace(0.2,1,length(iSat));

% Visualize melt fraction at each timestep
for i = 1:tn;

    % Extract melt fraction for the time step
    zr = zrAll(:,:,:,i); 
    
    % Time into the simulation to which this timestep corresponds
    t(i) = tstart+i*twriteout*dt;
    t(i) = t(i)/1e3; % [yr --> kyr]
    
    % Begin the 3D figure
    figure('Renderer','opengl')
    ah = gca;
    
    % Plot the 3D image countours
    for j=1:length(iSat)
        p(j) = patch(isosurface(X,Y,Z,zr(:,:,:),iSat(j)));
        isonormals(X,Y,Z,zr(:,:,:),p(j))
        set(p(j),'CData',j,'FaceAlpha',iTran(j));
    end
    
    % Tweak the countour color properties
    set(p(1:length(iSat)),'CDataMapping','direct','FaceColor','flat','EdgeColor','none')
    
    % Set the colormap
    cmap = flip(hot(length(iSat)));
    colormap(cmap)
    
    % Make a colobar and tweak it so iso-levels correspond to their color
    caxis([0 length(iSat)])
    cbh = colorbar('YTick',(1:length(iSat))+0.5, 'YTickLabel',num2str(iSat(:)));
    ylabel(cbh,'\DeltaZr (ppm)','FontSize',12);
    
    % Plot accoutrements
    axes(ah)
    axis equal
    
    % Define axis extents and label, keeping in mind the 3D orientation
    xlim([0 xdim]); ylim([0 zdim]); zlim([0 ydim])
    xlabel('distance X (km)','FontSize',12);
    ylabel('distance Z (km)','FontSize',12);
    zlabel('depth Y (km)','FontSize',12);
    % Give a title
    tstr = sprintf('Zircon Saturation\n%d ky',t(i));
    th = title(tstr,'FontSize',14);
    
    % Tweak the viewing geometry
    box on; grid on; daspect([1 1 1])
    view([-50,16]);
    camproj perspective
    camlight; lighting gouraud;
    rotate3d on

    % Save the figure as a tiff and MATLAB figure
    fname = sprintf('zrSaturation_F%d_R%d_T%02d',flux,reali,i);
    print(gcf,'-dtiffn',fname);
    %savefig(gcf,fname);
    close(gcf)
end