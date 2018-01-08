l=linspace(-1e-3,3e-3,100);
F=zeros(size(l));
for i=1:100
    if l(i)<0
        F(i)=0.3;
    else
        F(i)=-0.3*1e3*l(i)+0.3;
    end
    if F(i)<0
        F(i)=0;
    end
end
figure
plot(l,F)
axis([-1e-3 3e-3 0 0.4])
title('Pushing force vs. elongation model')
xlabel('$\Delta l [\mu m/s]$','interpreter','latex')
ylabel('Force [N]')
grid on