function bsc_plotClassifiedStreamsWithTube(wbFG, classification, t1, figView, saveDir, subSelect, colors, tubeSelect, tubeColors, tubeRadius)
%
%   bsc_plotClassifiedStreamsWithTube(wbFG, classification, t1, figView, saveDir, subSelect, colors, tubeSelect, tubeColors, tubeRadius)
%
%   PURPOSE: Same as bsc_plotClassifiedStreamsAdaptive_v2, but overlays a
%   centroid "tube" on top of the streamlines for the requested tracts.
%   Matches the MDLF-style figure: thin streamlines for context, thick
%   solid tube running along the tract core.
%
%  -wbFG, classification, t1, figView, saveDir, subSelect, colors:
%   See bsc_plotClassifiedStreamsAdaptive_v2 — identical semantics.
%
%  -tubeSelect: vector of tract indices (into classification.names) that
%   should also get a centroid tube. If empty or undefined, all tracts in
%   subSelect get a tube.
%
%  -tubeColors: Nx3 RGB matrix (rows aligned with tubeSelect). If
%   undefined, defaults to the corresponding streamline color.
%
%  -tubeRadius: scalar tube radius in mm. Defaults to 2.
%
% (C) Adapted from bsc_plotClassifiedStreamsAdaptive_v2 (Daniel Bullock, 2017).


%% preliminaries
[wbFG, fe] = bsc_LoadAndParseFiberStructure(wbFG);

if ischar(classification)
    load(classification);
end

if ~isempty(fe)
    if length(fe.life.fit.weights)==length(classification.index)
        classification = wma_clearNonvalidClassifications(classification, fe);
    else
        warning('mismatch between classification structure and fe weights')
    end
end

if ischar(t1)
    t1 = niftiRead(t1);
end

if notDefined('subSelect')
    subSelect = 1:length(classification.names);
end

if notDefined('colors')
    colors = zeros(length(subSelect), 3);
    if rem(length(subSelect),2)==0 && sum(subSelect(2:2:end)==subSelect(1:2:end)+1)==length(subSelect)/2
        for iTracts = 2:2:length(subSelect)
            c = rand(1,3);
            colors(iTracts,:)   = c;
            colors(iTracts-1,:) = c;
        end
    else
        for iTracts = 1:length(subSelect)
            colors(iTracts,:) = rand(1,3);
        end
    end
end

if notDefined('saveDir'),    saveDir = pwd;          end
if notDefined('tubeSelect'), tubeSelect = subSelect; end
if notDefined('tubeRadius'), tubeRadius = 2;         end

if notDefined('tubeColors')
    tubeColors = zeros(length(tubeSelect), 3);
    for k = 1:length(tubeSelect)
        match = find(subSelect == tubeSelect(k), 1);
        if isempty(match)
            tubeColors(k,:) = [1 1 1];
        else
            tubeColors(k,:) = colors(match,:);
        end
    end
end

%% plotting preliminaries
tractStruc = bsc_makeFGsFromClassification(classification, wbFG);

if length(subSelect)==1
    meanVec = zeros(length(tractStruc(subSelect).fg.fibers),3);
    for istreams = 1:length(tractStruc(subSelect).fg.fibers)
        meanVec(istreams,:) = mean(tractStruc(subSelect).fg.fibers{istreams},2);
    end
    totalMeans = mean(meanVec,1);
    slices = {[-1 0 0],[0 totalMeans(2) 0],[0 0 totalMeans(3)-10]};
else
    slices = {[-1 0 0],[0 -1 0],[0 0 -1]};
end

fh = figure('name','FiberClassificationPlotWithTube','color','k','units','normalized','position',[.5 .5 .5 .5]);
axis equal
fhNum = fh.Number;
hold on

switch lower(figView)
    case 'saggital'
        h = mbaDisplayBrainSlice(t1, slices{1});
    case 'coronal'
        h = mbaDisplayBrainSlice(t1, slices{2});
    case {'axial','transverse'}
        h = mbaDisplayBrainSlice(t1, slices{3});
end

%% streamlines
for itract = 1:length(subSelect)
    tractIndex = subSelect(itract);
    if exist('lh','var'), delete(lh); end
    if ~isempty(tractStruc(tractIndex).fg.fibers)
        [fh, lh] = mbaDisplayConnectome(tractStruc(tractIndex).fg.fibers, fhNum, colors(itract,:), 'single');
        delete(lh)
        fprintf('\n %s \n', classification.names{tractIndex})
    end
end

%% tube overlay
for k = 1:length(tubeSelect)
    tractIndex = tubeSelect(k);
    fg = tractStruc(tractIndex).fg;
    if isempty(fg.fibers), continue; end
    centroid = localComputeCentroid(fg, 100);
    if isempty(centroid), continue; end
    localDrawTube(centroid, tubeRadius, tubeColors(k,:), 24);
    fprintf('\n tube: %s \n', classification.names{tractIndex})
end

%% views
switch lower(figView)
    case 'saggital'
        fig.views   = {[90,0],[-90,0]};
        light.angle = {[45,15],[-45,15]};
        fig.names   = {strcat(wbFG.name,'_ClassifiedTracksTube_Saggital_',num2str(subSelect))};
    case 'coronal'
        fig.views   = {[0,0],[180,0]};
        light.angle = {[25,25],[155,-25]};
        fig.names   = {strcat(wbFG.name,'_ClassifiedTracksTube_Coronal_',num2str(subSelect))};
    case {'axial','transverse'}
        fig.views   = {[0,90],[0,-90]};
        light.angle = {[1,89],[1,-89]};
        fig.names   = {strcat(wbFG.name,'_ClassifiedTracksTube_Axial_',num2str(subSelect))};
end

for iview = 1:length(fig.views)
    fprintf('%s %s', figView, num2str(iview));
    view(fig.views{iview})

    if ischar(light.angle{iview})
        camlight(light.angle{iview})
    else
        lh = camlight(light.angle{iview}(1), light.angle{iview}(2));
    end

    figureHand = gcf;
    figureHand.Units = 'inches';

    fprintf('Saving as %s', strcat(fig.names{1},num2str(iview)))
    saveas(fh, strrep(strcat(saveDir,'/',fig.names{1},num2str(iview)),'  ','_'), 'png')
end

close all

end


function c = localComputeCentroid(fg, nPoints)
nF = length(fg.fibers);
R  = zeros(nF, nPoints, 3);
keep = false(nF,1);
for i = 1:nF
    f = fg.fibers{i}';
    d = [0; cumsum(sqrt(sum(diff(f).^2, 2)))];
    if d(end) < eps, continue; end
    s = linspace(0, d(end), nPoints);
    for k = 1:3
        R(i,:,k) = interp1(d, f(:,k), s, 'linear');
    end
    keep(i) = true;
end
R = R(keep,:,:);
if isempty(R), c = []; return; end

ref = squeeze(R(1,:,:));
for i = 2:size(R,1)
    cur   = squeeze(R(i,:,:));
    dFwd  = sum(sqrt(sum((cur          - ref).^2, 2)));
    dRev  = sum(sqrt(sum((flipud(cur)  - ref).^2, 2)));
    if dRev < dFwd
        R(i,:,:) = flipud(cur);
    end
end

c = squeeze(mean(R, 1));
end


function localDrawTube(curve, radius, color, nCircle)
N = size(curve,1);
T = diff(curve);
T = [T; T(end,:)];
T = T ./ max(sqrt(sum(T.^2,2)), eps);

seed = [0 0 1];
if abs(T(1,:) * seed') > 0.99, seed = [0 1 0]; end
n0 = cross(T(1,:), seed); n0 = n0 / norm(n0);

Nrm = zeros(N,3);
Nrm(1,:) = n0;
for i = 2:N
    a = T(i-1,:); b = T(i,:);
    ax = cross(a, b);
    s  = norm(ax);
    cs = dot(a, b);
    v  = Nrm(i-1,:);
    if s < 1e-8
        Nrm(i,:) = v;
    else
        ax  = ax / s;
        ang = atan2(s, cs);
        Nrm(i,:) = v*cos(ang) + cross(ax, v)*sin(ang) + ax*(ax*v')*(1-cos(ang));
    end
    Nrm(i,:) = Nrm(i,:) / norm(Nrm(i,:));
end

B = cross(T, Nrm, 2);
B = B ./ max(sqrt(sum(B.^2,2)), eps);

theta = linspace(0, 2*pi, nCircle);
X = zeros(N, nCircle); Y = X; Z = X;
for i = 1:N
    off = radius*(cos(theta)'*Nrm(i,:) + sin(theta)'*B(i,:));
    X(i,:) = curve(i,1) + off(:,1)';
    Y(i,:) = curve(i,2) + off(:,2)';
    Z(i,:) = curve(i,3) + off(:,3)';
end

surf(X, Y, Z, 'FaceColor', color, 'EdgeColor', 'none', ...
    'FaceLighting', 'gouraud', 'AmbientStrength', 0.4, ...
    'DiffuseStrength', 0.7, 'SpecularStrength', 0.3);
end
