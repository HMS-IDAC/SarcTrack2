function J = imDrawSarcomeresCB(I,rs1,cs1,as1,ds1,ds0) % CB: color bar

C = cool(size(I,1));
CB = zeros(size(I,1),20,3);
for i = 1:size(I,1)
    c = C(size(I,1)-i+1,:);
    CB(i,:,:) = repmat(reshape(c,[1 1 3]),[1 size(CB,2)]);
end
CB = [0.5*ones(size(I,1),5,3) CB];
C = cool(length(ds0));

rods0 = cell(1,length(rs1));
rods1 = cell(1,length(rs1));
rods2 = cell(1,length(rs1));
idss = zeros(1,length(rs1));
for i = 1:length(rs1)
    idss(i) = find(ds1(i) == ds0);
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
J = I;
J = insertShape(J,'Line',rods0,'Color','black');

for i = 1:length(ds0)
    idx = find(idss == i);
    if ~isempty(idx)
        J = insertShape(J,'Line',rods1(idx),'Color',C(i,:));
        J = insertShape(J,'Line',rods2(idx),'Color',C(i,:));
    end
end

fs = 8;
J = insertText(J,[cs1' rs1']-fs,cellstr(int2str((1:length(rs1))')),'BoxOpacity',0,'TextColor','white','FontSize',fs);
J = [J CB];
    
end