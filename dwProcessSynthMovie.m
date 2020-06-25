%% parameters

ds = 9:0.1:11;
stretch = 1; % stretch of morlet wavelet
scale = 1.5; % scale of morlet wavelet
nangs = 8;

%% read frames (from memory)

nFrames = size(S,3);

%% read frames (from disk)

% disp('reading frames')
% 
% v = VideoReader('~/Downloads/SarcTrackSampleVideos/Synth/Sample_04_12.avi');
% nFrames = round(v.Duration*v.FrameRate);
% S = zeros(v.Height,v.Width,nFrames);
% count = 0;
% while hasFrame(v)
%     if mod(count+1,round(nFrames/10)) == 1
%         fprintf('.')
%     end
%     count = count+1;
%     frame = readFrame(v);
%     I = double(rgb2gray(frame))/255;
%     S(:,:,count) = I;
% end
% fprintf('\n')
% tlvt(S)

%% initialize tracks

I = normalize(S(:,:,1));
[rs,cs,as,sp,~,K,imDA] = imFindSarcomeres(I,ds,nangs,stretch,scale);
imshowpair(I,K)

%% compute ridge-evidence volume

V = zeros(size(S));
cimDA = cell(1,size(S,3));

V(:,:,1) = K;
cimDA{1} = imDA;
for iFrame = 2:nFrames
    disp(iFrame)
    I = normalize(S(:,:,iFrame));

    [~,~,~,~,~,V(:,:,iFrame),cimDA{iFrame}] = imFindSarcomeres(I,ds,nangs,stretch,scale);
end

%% check detection

for i = 1:5:nFrames
    bw = imregionalmax(V(:,:,i).*imbinarize(V(:,:,i)));
    imDA = cimDA{i};
    [rs1,cs1] = find(bw);
    rs1 = rs1';
    cs1 = cs1';
    as1 = zeros(1,length(rs1));
    ds1 = zeros(1,length(rs1));
    for j = 1:length(rs1)
        ds1(j) = ds(imDA(rs1(j),cs1(j),1));
        as1(j) = (imDA(rs1(j),cs1(j),2)-1)/nangs*pi;
    end
    J = repmat(normalize(S(:,:,i)),[1 1 3]);
    J = imDrawSarcomeresCB(J,rs1,cs1,as1,ds1,ds);
    imshow(J)
    pause(0.1)
end

%% track

d = 7;
x = rs;
y = cs;
tracks = zeros(2,nFrames,length(x));
for i = 1:length(x)
    disp(i/length(x))

    row = x(i);
    col = y(i);

%     imshow(V(:,:,1)), hold on
%     plot(col,row,'o'), hold off

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
%         volumeViewer(P)
%         return

    else
        tracks(:,:,i) = repmat([row; col],[1 nFrames]);
    end
end

%% visualize tracks with angles, separations; gather data

nTracks = size(tracks,3);
rcasm = zeros(4,nFrames,nTracks); % rows, cols, angles, separations, magnitudes
for iFrame = 1:nFrames
    disp(iFrame)
    I = normalize(S(:,:,iFrame));
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
    J = imDrawSarcomeresCB(repmat(I,[1 1 3]),rs1,cs1,as1,ds1,ds);
    imshow(J)
    pause(0.1)
end

%% plot separations

for iTrack = 1:nTracks
    sp = rcasm(4,:,iTrack);
    plot(sp)
    axis([1 nFrames min(ds) max(ds)])
    pause(0.05)
end

%% plot average separations

plot(mean(rcasm(4,:,:),3)), axis([1 nFrames min(ds) max(ds)])

%% filter tracks based on proximity

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

imshow(mean(V,3)), hold on
plot(meanC,meanR,'.k')
plot(meanC(idx2remP),meanR(idx2remP),'xr'), hold off

%% filter tracks based on intensity

idx2remI = find(meanM < prctile(meanM,20))';
imshow(mean(V,3)), hold on
plot(meanC,meanR,'ok')
plot(meanC(idx2remI),meanR(idx2remI),'or'), hold off

%% filter tracks based on position variance

nTracks = size(tracks,3);
circs = [];
for iTrack = 1:nTracks
    disp(iTrack/nTracks)
    dr = diff(rcasm(1,:,iTrack));
    dc = diff(rcasm(2,:,iTrack));
    devs = sqrt(sum([dr.^2; dc.^2]));
    circs = [circs; [meanC(iTrack) meanR(iTrack) max(devs)]];
end
J = repmat(mean(V,3),[1 1 3]);
rads = circs(:,3);
thr = 3.5;
J = insertShape(J,'Circle',circs(rads <= thr,:),'Color','green');
J = insertShape(J,'Circle',circs(rads > thr,:),'Color','red');
imshow(J)
idx2remV = find(rads > thr)';

%% filter tracks based on angle variance (apply only to remaining tracks to save time)

% figure
idx2remA = [];
idx2remPIV = unique([idx2remP idx2remI idx2remV]);
idx2keep = 1:nTracks;
idx2keep(idx2remPIV) = [];
for iSelTrack = 1:length(idx2keep)
    iTrack = idx2keep(iSelTrack);
    disp(iTrack/nTracks)
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
imshow(mean(V,3)), hold on
plot(meanC,meanR,'ok')
plot(meanC(idx2remA),meanR(idx2remA),'or'), hold off


%% aggregate ind2rem

idx2rem = unique([idx2remP idx2remI idx2remV idx2remA]);
imshow(mean(V,3)), hold on
plot(meanC,meanR,'ok')
plot(meanC(idx2rem),meanR(idx2rem),'or'), hold off

%% visualize selected tracks with angles, separations

nTracks = size(tracks,3);
idx2keep = 1:nTracks;
idx2keep(idx2rem) = [];
stracksK = tracks(:,:,idx2keep);
stracksR = tracks(:,:,idx2rem);

outFit = '~/Downloads/Outputs';
if ~exist(outFit,'dir')
    mkdir(outFit);
end

for iFrame = 1:nFrames
    disp(iFrame)
    I = normalize(S(:,:,iFrame));
    imDA = cimDA{iFrame};
    M = V(:,:,iFrame);
    
    rs1 = squeeze(stracksK(1,iFrame,:))';
    cs1 = squeeze(stracksK(2,iFrame,:))';
    as1 = zeros(1,length(rs1));
    ds1 = zeros(1,length(rs1));
    for j = 1:length(rs1)
        ds1(j) = ds(imDA(rs1(j),cs1(j),1));
        as1(j) = (imDA(rs1(j),cs1(j),2)-1)/nangs*pi;
    end
    J = imDrawSarcomeresCB(repmat(I,[1 1 3]),rs1,cs1,as1,ds1,ds);
    
%     rs1 = squeeze(stracksR(1,iFrame,:))';
%     cs1 = squeeze(stracksR(2,iFrame,:))';
%     as1 = zeros(1,length(rs1));
%     ds1 = zeros(1,length(rs1));
%     for j = 1:length(rs1)
%         ds1(j) = ds(imDA(rs1(j),cs1(j),1));
%         as1(j) = (imDA(rs1(j),cs1(j),2)-1)/nangs*pi;
%     end
%     J = imDrawSarcomeres(J,rs1,cs1,as1,ds1,'red');
    
    imshow(J)
    pause(0.1)
    imwrite(J,[outFit filesep sprintf('Frame%05d.tif',iFrame)]);
end

%% plot average displacements for selected tracks

plot(mean(rcasm(4,:,idx2keep),3),'g'), hold on
plot(mean(rcasm(4,:,:),3),'k'), hold off, axis([1 nFrames min(ds) max(ds)])


%% estimate frequency
disp('estimating frequency')

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

figureQSS
plot(fpd.x,fpd.y2Fit,'r'), hold on
plot(fpd.x,fpd.ySin,'g'),
plot(fpd.x,fpd.ySaw,'b'),
plot(fpd.x,fpd.ySawFit,'k'), hold off
legend('y2Fit','ySin','ySaw', 'ySawFit')


%% fit sawtooth curves
disp('fitting sawtooth curves')

prms = [];
idcs = [];
dsls = [];
figureQSS
for index = 1:length(idx2keep)
    if mod(index,round(length(idx2keep)/10)) == 1
        fprintf('.')
    end

    x = (0:nFrames-1)';
    dsl = rcasm(4,:,idx2keep(index))';
    
%     plot(x,dsl)
%     pause(0.1)
%     continue
    
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
            idcs = [idcs index];
            dsls = [dsls dsl];
            
            plot(fpd.x,y,'r'), hold on
            plot(fpd.x,ySawFit,'g'), hold off
            legend('y2Fit','ySawFit')
            pause(0.1)
        end
    end
end
close all
fprintf('\n')

%% plot stats

prms0 = prms;
if ~isempty(prms0)
    prd = 2*pi/f.b1;
    c = prms0(1,:)/(2*pi)*prd;
    r = prms0(2,:)/(2*pi)*prd;
    o = prms0(3,:)/(2*pi)*prd;
    
    T = array2table([c' r' o' prms0(4:7,:)'],'VariableNames',{'contraction_time','relaxation_time','offset_from_average','min_ds','max_ds','min_ds_fit','max_ds_fit'});
    A = table2array(T);
    prms = A';
    labels = cell(1,7*size(prms,2));
    for i = 1:size(prms,2)
        labels{i} = 'contraction time';
        labels{size(prms,2)+i} = 'relaxation time';
        labels{2*size(prms,2)+i} = 'offset from average';
        labels{3*size(prms,2)+i} = 'min ds';
        labels{4*size(prms,2)+i} = 'max ds';
        labels{5*size(prms,2)+i} = 'min ds fit';
        labels{6*size(prms,2)+i} = 'max ds fit';
    end
    figureQSS
    boxplot([prms(1,:) prms(2,:) prms(3,:)],labels(1:3*size(prms,2)))
    figureQSS
    boxplot([prms(4,:) prms(5,:) prms(6,:) prms(7,:)],labels(3*size(prms,2)+1:end))

    % -------------------------
    % histograms

    figureQSS
    titles = {'contraction time','relaxation time','offset from average','min ds','max ds','min ds fit','max ds fit'};
    for i = 1:7
        subplot(1,7,i)
        histogram(prms(i,:))
        title(titles{i})
        fprintf('median %s: %f\n', titles{i}, median(prms(i,:)));
    end
    fprintf('median min/max: %f\n', median(prms(4,:)./prms(5,:)));
    fprintf('median min/max fit: %f\n', median(prms(6,:)./prms(7,:)));
    
    
    
    vn = cell(1,size(dsls,2));
    for i = 1:size(dsls,2)
        vn{i} = sprintf('track%05d',idcs(i));
    end
    T = array2table(dsls,'VariableNames',vn);
    A = table2array(T);
    figure
    histogram(A(:),20), title('dists')

    T = array2table([prd 1/prd],'VariableNames',{'period','frequency'});
    disp(T)

end

%% plot average fit curve, compare to ground truth

hPad = 10;
p = ones(hPad,1);
t = [p*A(1,1); A(:,1); p*A(end,1)];
at = t;
for i = 2:size(A,2)
    cs = zeros(1,2*hPad);
    for j = 1:2*hPad
%         plot(t,'r'), hold on
        u = [A(1,i)*ones(j,1); A(:,i); A(end,i)*ones(2*hPad-j,1)];
%         plot(u,'g'), hold off
%         pause
        cs(j) = corr(t,u);
    end
    [~,j] = max(cs);
    u = [A(1,i)*ones(j,1); A(:,i); A(end,i)*ones(2*hPad-j,1)];
%     plot(t,'r'), hold on
%     plot(u,'g'), hold off
%     pause
    at = at+u;
end
at = at/size(A,2);
% plot(at)

% uncomment the next two lines if reading movie from disk
% [~,vName] = fileparts(v.Name);
% load([v.Path filesep vName '.mat']); % loads 'dsts'

i0 = round((length(at)-length(dsts))/2);
at = at(i0:i0+length(dsts)-1);

t = [p*dsts(1); dsts'; p*dsts(end)];
cs = zeros(1,2*hPad);
for j = 1:2*hPad
    u = [at(1)*ones(j,1); at; at(end)*ones(2*hPad-j,1)];
    cs(j) = corr(t,u);
end
[~,j] = max(cs);
u = [at(1)*ones(j,1); at; at(end)*ones(2*hPad-j,1)];

plot(t,'g'), hold on
plot(u,'b'), hold off
axis([-1 length(t)+1 9.25 10.75])
% legend(sprintf('gt (%d, %d)', f1, f2),sprintf('avg fit (%d)',size(A,2)))
title(sprintf('avg fit from %d measurements',size(A,2)))