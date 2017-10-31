%% retrieve experimental data
adhesion_dataset_1
xdata=deltal;
ydata=adhesion_force;
%% data fitting
[x,sigma]=lsqcurvefit(@(x,xdata) x(1)*xdata.*exp(-x(2)*(xdata.^x(3))),ones(3,1),xdata,ydata);
%% plot
hold on
elongation=linspace(0,2,100);
plot(elongation,x(1)*elongation.*exp(-x(2)*(elongation.^x(3))), 'linewidth',2)
hold off