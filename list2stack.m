function S = list2stack(l)

S = zeros(size(l{1},1),size(l{1},2),length(l));
for i = 1:length(l)
    S(:,:,i) = l{i};
end

end