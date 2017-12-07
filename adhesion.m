clearvars
load 'adhesion'

X_adh     = zeros(3,3);     %this is a matrix of 0s, only for initialization purpus
sigma_adh = zeros(3,3);     %this is a matrix of 0s, only for initialization purpus

X     = zeros(5,3);         %this is a matrix of 0s, only for initialization purpus
sigma = zeros(5,1);         %this is a matrix of 0s, only for initialization purpus

for n = 1:3                     %this is a for cycle
    for m = 1:5                 %this is anather for cycl, wich is nested
        [X(m,:), sigma(m)] = lsqcurvefit(@(x,xdata) x(1)*xdata...
            .*exp(-x(2)*(xdata.^x(3))),ones(3,1),max(set(n).exp(m).el,zeros(size(set(n)...
            .exp(m).el))),max(set(n).exp(m).F,zeros(size(set(n).exp(m).F))));
    end
    X_adh(n,:) = sum(X)/5;
    sigma_adh(n,:) = sqrt(sum((X-X_adh(n,:)).^2)/4);
end
