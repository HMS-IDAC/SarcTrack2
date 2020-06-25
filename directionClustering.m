function [c,l] = directionClustering(V,k,clusterProximityThreshold,ignoreAngleSign)

% [c,l] = directionClustering(V,k,clusterProximityThreshold,ignoreAngleSign)
% clusters directions (angles), possibly ignoring sign (i.e. outputing clusters
% of angles in [0,pi) only)
%
% ------
% inputs
% ------
%
% V: nx2 matrix of points *in the circle*
% k: can be thought of as inverse of 'bandwidth' in mean shift literature parlance
%    rule of thumb: for n equally spaced clusters on the 2D disk set k = n
% clusterProximityThreshold: if the algorithm finds two cluster centers c1 and c2
%    and dot(c1,c2) > clusterProximityThreshold, one of the centers will not be considered
% ignoreAngleSign: true or false; if true, cluster centers will be in the range [0,pi),
%    otherwise [0,2*pi); use ignoreAngleSign = true when clustering 'lines',
%    and false when clustering 'arrows'.
%
% -------
% outputs
% -------
% c: cluster centers, one per row
% l: vector with indices of clusters;
%    e.g.: l(i) = j means point V(i,:) belongs to cluster c(j,:)
%
% -------
% example
% -------
%
% a = pi/4;
% n = 20;
% std = pi/12;
% A = a+std*randn(n,1);
% V = [cos(A) sin(A)];
% 
% a = a+pi;
% n = 20;
% std = pi/12;
% A = a+std*randn(n,1);
% V = [V; [cos(A) sin(A)]];
% 
% k = 2;
% clusterProximityThreshold = 0.9;
% ignoreAngleSign = true;
% [c,l] = directionClustering(V,k,clusterProximityThreshold,ignoreAngleSign);
% 
% figure, hold on
% for j = 1:size(c,1)
%     h = (j-1)/size(c,1)*2/3;
%     plot(V(l == j,1),V(l == j,2),'.','Color',hsv2rgb([h 1 1]))
%     plot(c(j,1),c(j,2),'s','Color',hsv2rgb([h 1 1]))
%     plot(c(j,1),c(j,2),'.k')
% end
% axis equal, axis([-1.1 1.1 -1.1 1.1]), grid on
% hold off
% title(sprintf('%d clusters found', size(c,1)))
%
% ---------
% reference
% ---------
% T. Kobayashi and N. Otsu, "Von Mises-Fisher Mean Shift for Clustering on a Hypersphere,"
% Pattern Recognition (ICPR), 2010 20th International Conference on, Istanbul, 2010, pp. 2130-2133.
%
% -----------
% development
% -----------
% Marcelo Cicconet, Sep 2017

if ignoreAngleSign
%     plot(V(:,1),V(:,2),'.'), hold on
%     axis equal, axis([-1.1 1.1 -1.1 1.1]), grid on

    as = atan2(V(:,2),V(:,1)); % [-pi, pi]
    as(as < 0) = 2*pi+as(as < 0); % [0, 2*pi]

    as(as >= pi) = as(as >= pi)-pi; % when angle > pi direction is equivalent to angle-pi

%     W = [cos(as) sin(as)];
%     figure
%     plot(W(:,1),W(:,2),'o'), hold off
%     axis equal, axis([-1.1 1.1 -1.1 1.1]), grid on

    as = 2*as; % [0,pi] -> [0,2*pi] for mean shift 'continuity'
    W = [cos(as) sin(as)];

%     figure
%     plot(W(:,1),W(:,2),'o'), hold off
%     axis equal, axis([-1.1 1.1 -1.1 1.1]), grid on
else
    W = V;
end
    
% mean shift

[n,d] = size(W);
K = @(x,x0) exp(k*dot(x,x0)); % kernel

C = zeros(n,d); % convergence points for every sample
for i = 1:n
    m = W(i,:);

    previousm = -m;
    while dot(previousm,m) < 0.999999
        previousm = m;

        M = zeros(1,d);

        nbh = find(sum(W.*repmat(m,[n 1]),2) > 0.7);

        for jj = 1:length(nbh)
            j = nbh(jj);
            M = M+K(W(j,:),m)*W(j,:);
        end

        m = M/norm(M);
    end

    C(i,:) = m;
end

c = C(1,:); % convergence points (one per cluster)
l = zeros(n,1); % cluster labels (one per sample)
for i = 1:n
    didbreak = 0;
    for j = 1:size(c,1)
        if dot(C(i,:),c(j,:)) > clusterProximityThreshold
            l(i) = j;
            didbreak = 1;
            break;
        end
    end
    if ~didbreak
        c = [c; C(i,:)];
        l(i) = size(c,1);
    end
end

if ignoreAngleSign
    as = atan2(c(:,2),c(:,1)); % [-pi, pi]
    as(as < 0) = 2*pi+as(as < 0); % [0, 2*pi]
    as = as/2; % back to [0, 180]
    c = [cos(as) sin(as)];
end

end