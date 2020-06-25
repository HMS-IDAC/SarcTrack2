%% parameters

frameIndex = 53;

ds = 9:0.1:11; % range of distances (in pixels)
stretch = 1; % stretch of morlet wavelet
scale = 1.5; % scale of morlet wavelet
nangs = 8;

%% read frame, show wavelets

I = S(:,:,frameIndex);

mr1 = normalize(smorlet2(ds(1),stretch,scale,90));
mr2 = normalize(smorlet2(ds(end),stretch,scale,90));
J = normalize(I);
J(1:size(mr1,1),1:size(mr1,2)) = mr1;
J(size(mr1,1)+1:size(mr1,1)+size(mr2,1),1:size(mr2,2)) = mr2;
imtool(J)

%% process

[rs,cs,as,sp,~,K,~] = imFindSarcomeres(I,ds,nangs,stretch,scale);
imshowpair(I,K)

%% draw

J = imDrawSarcomeresCB(repmat(I,[1 1 3]),rs,cs,as,sp,ds);
imshow(J)