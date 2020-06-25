function dwProcessVideo(path)

% clear, clc
tic

%% parameters

ds = 9:0.1:11;
stretch = 1; % stretch of morlet wavelet
scale = 1.5; % scale of morlet wavelet
nangs = 8; % number of rotation angles

%% read frames
disp('reading frames')

v = VideoReader(path);
nFrames = round(v.Duration*v.FrameRate);
S = zeros(v.Height,v.Width,nFrames);
count = 0;
while hasFrame(v)
    if mod(count+1,round(nFrames/10)) == 1
        fprintf('.')
    end
    count = count+1;
    frame = readFrame(v);
    I = double(rgb2gray(frame))/255;
    S(:,:,count) = I;
end
fprintf('\n')


%% initialize tracks
disp('initialize tracks')

I = normalize(S(:,:,1));
[rs,cs,~,~,~,K,imDA,W] = imFindSarcomeres(I,ds,nangs,stretch,scale);

%% compute ridge-evidence volume
disp('compute ridge-evidence volume')

V = zeros(size(S));
cimDA = cell(1,size(S,3));

V(:,:,1) = K;
cimDA{1} = imDA;
fprintf('.')
for iFrame = 2:nFrames
    if mod(iFrame,round(nFrames/10)) == 1
        fprintf('.')
    end
    I = normalize(S(:,:,iFrame));
    [~,~,~,~,~,V(:,:,iFrame),cimDA{iFrame}] = imFindSarcomeres(I,ds,nangs,stretch,scale,W);
end
fprintf('\n')

%% track
disp('tracking via dynamic programming')

d = 7;
x = rs;
y = cs;
tracks = zeros(2,nFrames,length(x));
for i = 1:length(x)
    if mod(i,round(length(x)/10)) == 1
        fprintf('.')
    end

    row = x(i);
    col = y(i);

%     imshow(V(:,:,1)), hold on
%     plot(col,row,'o'), hold off
%     return

    if row-d >= 1 && row+d <= size(I,1) && col-d >= 1 && col+d <= size(I,2)
        CV = V(row-d:row+d,col-d:col+d,:);
    
        CV = CV-min(CV(:));
        CV = CV/max(CV(:));
        C = 1-CV;

        ijPath = dpV(C,d+1,d+1);
        tracks(:,:,i) = ijPath-(d+1)*ones(size(ijPath))+repmat([row; col],[1 size(ijPath,2)]);

%         P = 0.5*CV;
%         for k = 1:size(ijPath,2)
%             P(ijPath(1,k),ijPath(2,k),k) = 1;
%         end
%         tlvt(P)
%         return

    else
        tracks(:,:,i) = repmat([row; col],[1 nFrames]);
    end
end
fprintf('\n')

%% gather data
disp('gather data')

nTracks = size(tracks,3);
rcasm = zeros(4,nFrames,nTracks); % rows, cols, angles, separations, magnitudes
for iFrame = 1:nFrames
    imDA = cimDA{iFrame};
    M = V(:,:,iFrame);
    
    rs1 = squeeze(tracks(1,iFrame,:))';
    cs1 = squeeze(tracks(2,iFrame,:))';
    as1 = zeros(1,length(rs1));
    ds1 = zeros(1,length(rs1));
    ms1 = zeros(1,length(rs1));
    for j = 1:length(rs1)
        ds1(j) = ds(imDA(rs1(j),cs1(j),1));
        as1(j) = (imDA(rs1(j),cs1(j),2)-1)/nangs*pi;
        ms1(j) = M(rs1(j),cs1(j));
    end
    rcasm(1,iFrame,:) = rs1;
    rcasm(2,iFrame,:) = cs1;
    rcasm(3,iFrame,:) = as1;
    rcasm(4,iFrame,:) = ds1;
    rcasm(5,iFrame,:) = ms1;
end

%% filter tracks based on proximity
disp('filter tracks based on proximity')

meanR = squeeze(mean(rcasm(1,:,:),2));
meanC = squeeze(mean(rcasm(2,:,:),2));
meanM = squeeze(mean(rcasm(5,:,:),2));

dist = squareform(pdist([meanR,meanC]));
dist(triu(ones(size(dist))) > 0) = Inf;

idx2remP = [];

while 1
    [ii,jj] = find(dist == min(dist(dist~=0)));
    i = ii(1); j = jj(1);

    if dist(i,j) < 0.5*mean(ds)
        if meanM(i) < meanM(j)
            idx2remP = [idx2remP i];
        else
            idx2remP = [idx2remP j];
        end
        dist(i,j) = Inf;
        
%         plot(meanR,meanC,'.g'), hold on
%         plot(meanR([i j]),meanC([i j]),'.r'), hold off
%         pause
    else
        break
    end
end

% imshow(mean(V,3)), hold on
% plot(meanC,meanR,'.k')
% plot(meanC(idx2remP),meanR(idx2remP),'xr'), hold off

%% filter tracks based on intensity
disp('filter tracks based on intensity')

idx2remI = find(meanM < prctile(meanM,20))';
% imshow(mean(V,3)), hold on
% plot(meanC,meanR,'ok')
% plot(meanC(idx2remI),meanR(idx2remI),'or'), hold off

%% filter tracks based on position variance
disp('filter tracks based on position variance')

nTracks = size(tracks,3);
circs = [];
for iTrack = 1:nTracks
    dr = diff(rcasm(1,:,iTrack));
    dc = diff(rcasm(2,:,iTrack));
    devs = sqrt(sum([dr.^2; dc.^2]));
    circs = [circs; [meanC(iTrack) meanR(iTrack) max(devs)]];
end
rads = circs(:,3);
thr = 4;%3.5;
idx2remV = find(rads > thr)';
% J = repmat(mean(V,3),[1 1 3]);
% J = insertShape(J,'Circle',circs(rads <= thr,:),'Color','green');
% J = insertShape(J,'Circle',circs(rads > thr,:),'Color','red');
% imshow(J)


%% filter tracks based on angle variance (apply only to remaining tracks to save time)
disp('filter tracks based on angle variance')

% figure
idx2remA = [];
idx2remPIV = unique([idx2remP idx2remI idx2remV]);
idx2keep = 1:nTracks;
idx2keep(idx2remPIV) = [];
for iSelTrack = 1:length(idx2keep)
    if mod(iSelTrack,round(length(idx2keep)/10)) == 1
        fprintf('.')
    end
    
    iTrack = idx2keep(iSelTrack);
    as = rcasm(3,:,iTrack);
    x = cos(as);
    y = sin(as);
    xy = [x',y'];

    k = 2;
    clusterProximityThreshold = 0.9;
    ignoreAngleSign = true;
    c = directionClustering(xy,k,clusterProximityThreshold,ignoreAngleSign);

    if size(c,1) > 1
        idx2remA = [idx2remA iTrack];
    end
end
fprintf('\n')
% imshow(mean(V,3)), hold on
% plot(meanC,meanR,'ok')
% plot(meanC(idx2remA),meanR(idx2remA),'or'), hold off


%% aggregate ind2rem

disp('tracks flagged for removal:')
fprintf('proximity: %d, intensity: %d, position variance: %d, angle variance: %d\n', ...
    length(idx2remP), length(idx2remI), length(idx2remV), length(idx2remA));
idx2rem = unique([idx2remP idx2remI idx2remV idx2remA]);
% imshow(mean(V,3)), hold on
% plot(meanC,meanR,'ok')
% plot(meanC(idx2rem),meanR(idx2rem),'or'), hold off

%% estimate frequency
disp('estimating frequency')

nTracks = size(tracks,3);
idx2keep = 1:nTracks;
idx2keep(idx2rem) = [];
avgDsts = mean(rcasm(4,:,idx2keep),3)';

x = (0:nFrames-1)';
mx = prctile(avgDsts,90);
mn = prctile(avgDsts,10);
yAvg = (avgDsts-mean(avgDsts))/std(avgDsts);
y2Fit = 2*((avgDsts-mn)/(mx-mn)-0.5);

syAvg = smooth(yAvg,15);
f = fit(x,syAvg,'sin1');
ySin = f.a1*sin(f.b1*x+f.c1);

c = pi/2;
r = pi/2;
o = 0;
ySaw = f.a1*sawtooth(f.b1*x+pi/2+f.c1,c,r);

f2m = @(cro) -corr( f.a1*sawtooth(f.b1*x+cro(3)+pi/2+f.c1,cro(1),cro(2)) , y2Fit );

cro0 = [c; r; o];
lb = [0 0 -pi/4];
ub = [pi pi pi/4];
cro = fmincon(f2m,cro0,[],[],[],[],lb,ub,[],optimoptions('fmincon','Display','off'));
ySawFit = f.a1*sawtooth(f.b1*x+cro(3)+pi/2+f.c1,cro(1),cro(2));

% f plot data, used by dwCheckResults
fpd.x = x;
fpd.y2Fit = y2Fit;
fpd.ySin = ySin;
fpd.ySaw = ySaw;
fpd.ySawFit = ySawFit;
[fpdPath,fpdName] = fileparts(path);
save([fpdPath filesep fpdName '_fpd.mat'], 'fpd');


%% fit sawtooth curves
disp('fitting sawtooth curves')

prms = [];
idcs = [];
dsls = [];
for index = 1:length(idx2keep)
    if mod(index,round(length(idx2keep)/10)) == 1
        fprintf('.')
    end
    
    x = (0:nFrames-1)';
    dsl = rcasm(4,:,idx2keep(index))';
    
    mx = prctile(dsl,90);
    mn = prctile(dsl,10);
    if mx > mn
        y = 2*((dsl-mn)/(mx-mn)-0.5);
        f2m = @(prm) -corr( f.a1*sawtooth(f.b1*x+prm(3)+pi/2+f.c1,prm(1),prm(2)) , y );
        options = optimoptions('fmincon','Display','off');
        prm = fmincon(f2m,cro,[],[],[],[],lb,ub,[],options);
        ySawFit = f.a1*sawtooth(f.b1*x+prm(3)+pi/2+f.c1,prm(1),prm(2));
        z = ((ySawFit/2)+0.5)*(mx-mn)+mn; % original range

        rmse = sqrt(sum((y-ySawFit).^2)/nFrames);

        if rmse < 1 && max(dsl) < max(ds) && min(dsl) > min(ds)
            prms = [prms [prm; min(dsl); max(dsl); min(z); max(z)]];
            idcs = [idcs idx2keep(index)];
            dsls = [dsls dsl];
        end
    end
end
fprintf('\n')

%% draw outputs
disp('writing images')

[rpath,fname] = fileparts(path);
outFit = [rpath filesep fname '_DWFit'];
if ~exist(outFit,'dir')
    mkdir(outFit);
end

stracks = tracks(:,:,idcs);

for iFrame = 1:nFrames
    if mod(iFrame,round(nFrames/10)) == 1
        fprintf('.')
    end
    
    I = normalize(S(:,:,iFrame));
    imDA = cimDA{iFrame};
    
    rs1 = squeeze(stracks(1,iFrame,:))';
    cs1 = squeeze(stracks(2,iFrame,:))';
    as1 = zeros(1,length(rs1));
    ds1 = zeros(1,length(rs1));
    for j = 1:length(rs1)
        ds1(j) = ds(imDA(rs1(j),cs1(j),1));
        as1(j) = (imDA(rs1(j),cs1(j),2)-1)/nangs*pi;
    end
    J = imDrawSarcomeresCB(repmat(0.5*I,[1 1 3]),rs1,cs1,as1,ds1,ds);
    imwrite(J,[outFit filesep sprintf('Frame%03d.png',iFrame)]);
end
fprintf('\n')

%% write tables
disp('writing tables')

[rpath,fname] = fileparts(path);
outPathS = [rpath filesep fname '_DWStats.csv'];
outPathD = [rpath filesep fname '_DWDists.csv'];
outPathPF = [rpath filesep fname '_DWPrdFrq.csv'];


if ~isempty(prms)
    prd = 2*pi/f.b1;
    
    c = prms(1,:)/(2*pi)*prd;
    r = prms(2,:)/(2*pi)*prd;
    o = prms(3,:)/(2*pi)*prd;

    T = array2table([c' r' o' prms(4:7,:)'],'VariableNames',{'contraction_time','relaxation_time','offset_from_average','min_ds','max_ds','min_ds_fit','max_ds_fit'});
    writetable(T,outPathS);
    
    vn = cell(1,size(dsls,2));
    for i = 1:size(dsls,2)
        vn{i} = sprintf('track%05d',i);
    end
    T = array2table(dsls,'VariableNames',vn);
    writetable(T,outPathD);
    
    T = array2table([prd 1/prd],'VariableNames',{'period','frequency'});
    writetable(T,outPathPF);
else
    writetable(array2table([]),outPathS);
    writetable(array2table([]),outPathD);
end

%%
toc

end