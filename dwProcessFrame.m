clear, clc

%% read and display

path = '~/Downloads/SarcTrackSampleVideos/Synth/Sample_10_12.avi';

v = VideoReader(path);
T = zeros(v.Height,v.Width,v.FrameRate*v.Duration);
frameCount = 0;
while hasFrame(v)
    frameCount = frameCount+1;
    frame = readFrame(v);
    I = double(rgb2gray(frame))/255;
    T(:,:,frameCount) = I;
    disp(frameCount)
end
timeLapseViewTool(T);

%% parameters

frameIndex = 1;

ds = 9:0.1:11; % range of distances (in pixels)
stretch = 1; % stretch of morlet wavelet
scale = 1.5; % scale of morlet wavelet
nangs = 8;

%% read frame, show wavelets

v.CurrentTime = frameIndex/v.FrameRate;
frame = readFrame(v);
I = double(rgb2gray(frame))/255;

mr1 = normalize(smorlet2(ds(1),stretch,scale,90));
mr2 = normalize(smorlet2(ds(end),stretch,scale,90));
J = normalize(I);
J(1:size(mr1,1),1:size(mr1,2)) = mr1;
J(size(mr1,1)+1:size(mr1,1)+size(mr2,1),1:size(mr2,2)) = mr2;
imtool(J)

%% process

[rs,cs,as,sp,J,K,imDA] = imFindSarcomeres(I,ds,nangs,stretch,scale);
figureQSS
subplot(1,2,1), imshowpair(I,K)
subplot(1,2,2), imshow(J)
fprintf('# sarcomeres detected: %d\n',length(rs))

%% draw

J = imDrawSarcomeresCB(repmat(I,[1 1 3]),rs,cs,as,sp,ds);
imshow(J)