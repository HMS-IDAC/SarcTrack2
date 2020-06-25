function [rs,cs,as,sp,J,K,imDA,W] = imFindSarcomeres(I,ds,nangs,stretch,scale,varargin)
% W = varargin{1} = {wf,wt,wlist}

orientations = (0:nangs-1)*180/nangs;

% double-wavelet filters and templates

% ds = 11:0.2:13;
ndists = length(ds);
if nargin > 5
    wf = varargin{1}{1};
    wt = varargin{1}{2};
else
    wf = cell(ndists,nangs); % wavelet filters
    wt = cell(ndists,nangs); % wavelet templates
    for i = 1:ndists
        for j = 1:nangs
            mr = smorlet2(ds(i),stretch,scale,orientations(j));
            mrt = mr.*(mr > 0);
            wf{i,j} = mr;
            wt{i,j} = mrt;
    %         mr = imresize(mr,10,'nearest');
    %         mrt = imresize(mrt,10,'nearest');
    %         imshow([mr mrt],[])
    %         pause(0.5)
        end
    end
end

% convolutions

M = zeros(size(I,1),size(I,2),ndists,nangs);

% for i = 1:ndists
%     disp(i/ndists)
%     for j = 1:nangs
%         M(:,:,i,j) = conv2(I,wf{i,j},'same');
%     end
% end

mlist = cell(1,ndists*nangs);
if nargin > 5
    wlist = varargin{1}{3};
else
    wlist = cell(1,ndists*nangs);
    c = 0;
    for i = 1:ndists
        for j = 1:nangs
            c = c+1;
            wlist{c} = gpuArray(wf{i,j});
        end
    end
end
gpuI = gpuArray(I);
parfor i = 1:ndists*nangs
    mlist{i} = gather(conv2(gpuI,wlist{i},'same'));
end
c = 0;
for i = 1:ndists
    for j = 1:nangs
        c = c+1;
        M(:,:,i,j) = mlist{c};
    end
end

% regional maxima

[ma,ima] = max(M,[],4);
[md,imd] = max(ma,[],3);
bw = imregionalmax(md).*imbinarize(md);
% switchBetween(md,bw)

% detection

[rs0,cs0] = find(bw);
% rods = cell(1,length(rs0));
idd = zeros(1,length(rs0));
ida = zeros(1,length(rs0));
for i = 1:length(rs0)
    r = rs0(i);
    c = cs0(i);
    idxoptdist = imd(r,c);
    idxoptangl = ima(r,c,idxoptdist);
    idd(i) = idxoptdist;
    ida(i) = idxoptangl;
    
%     angle = orientations(idxoptangl)/180*pi;
%     d = ds(idxoptdist)/2*[cos(angle) sin(angle)];
%     xy = [r c];
%     rods{i} = [xy(2)-d(2) xy(1)-d(1) xy(2)+d(2) xy(1)+d(1)];
end
% J0 = insertShape(repmat(I,[1 1 3]),'Line',rods);
% imshow(J0)

%

% selection

rs1 = []; % rows
cs1 = []; % cols
as1 = []; % angles
ds1 = []; % distances
for idst = 1:ndists
%     disp(idst/ndists)
    for iang = 1:nangs
        idx = find(idd == idst & ida == iang);
        
        r = rs0(idx);
        c = cs0(idx);

        hks = floor(size(wf{idst,iang},1)/2);
        t = wt{idst,iang};
        t = t(:);
        t = t/sum(t);

        for i = 1:length(r)
            r0 = r(i)-hks;
            c0 = c(i)-hks;
            r1 = r(i)+hks;
            c1 = c(i)+hks;
            if r0 >= 1 && r1 <= size(I,1) && c0 >= 1 && c1 <= size(I,2)
                P = I(r0:r1,c0:c1);
                P = P/sum(P(:));
                if corr(t,P(:)) > 0.25
                    rs1 = [rs1 r(i)];
                    cs1 = [cs1 c(i)];
                    as1 = [as1 orientations(iang)/180*pi];
                    ds1 = [ds1 ds(idst)];
                end
            end
        end
    end
end

% draw

if nargout > 4
    rods0 = cell(1,length(rs1));
    rods1 = cell(1,length(rs1));
    rods2 = cell(1,length(rs1));
    for i = 1:length(rs1)
        d = ds1(i)/2*[cos(as1(i)) sin(as1(i))];
        xy = [rs1(i) cs1(i)];
        p = [xy(2)-d(2) xy(1)-d(1)];
        q = [xy(2)+d(2) xy(1)+d(1)];
        rods0{i} = [p q];
        d = ds1(i)/2*[cos(as1(i)+pi/2) sin(as1(i)+pi/2)];
        xy = fliplr(p);
        p2 = [xy(2)-d(2) xy(1)-d(1)];
        q2 = [xy(2)+d(2) xy(1)+d(1)];
        rods1{i} = [p2 q2];
        xy = fliplr(q);
        p2 = [xy(2)-d(2) xy(1)-d(1)];
        q2 = [xy(2)+d(2) xy(1)+d(1)];
        rods2{i} = [p2 q2];
    end
    J = repmat(normalize(I),[1 1 3]);
    J = insertShape(J,'Line',rods0,'Color','yellow');
    J = insertShape(J,'Line',rods1,'Color','green');
    J = insertShape(J,'Line',rods2,'Color','green');
end

% output

rs = rs1; % rows
cs = cs1; % cols
as = as1; % angles
sp = ds1; % separation between z-discs
if nargout > 5
    K = normalize(md);
end
if nargout > 6
    imDA = zeros(size(I,1),size(I,2),2);
    imDA(:,:,1) = imd;
    for i = 1:size(I,1)
        for j = 1:size(I,2)
            imDA(i,j,2) = ima(i,j,imd(i,j));
        end
    end
end
if nargout > 7
    W = {wf,wt,wlist};
end

end