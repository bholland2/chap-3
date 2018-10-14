function respiration1
% Parameter values
close all
v=1; % Vmax for use of substrate 
k=400; % Value of substrate at Vmax/2 
cA=20; % Initial A
cx0=0.2; % Shading rate 
nA=20; 
nx0=0.2; 
r1=0;
r2=0;

y1_0=0.0001;
y2_0=0.0001;
y3_0=0.1;
y4_0=0.1;
tspan=(0:1:1000);
%solver 
[T,Y]=ode15s(@thorn,tspan,[y1_0 y2_0 y3_0 y4_0],'RelTol'); % solver
carbon=Y(:,1);
nitrogen=Y(:,2);
leaf=Y(:,3);
root=Y(:,4);
SR=leaf./root;
plant=leaf+root;
RGR=gradient(log(plant));
kc=f(leaf);
kn=g(root);
rgr1=h(carbon,nitrogen);
rgr2=d(nitrogen,carbon);
g1=p(carbon,nitrogen,leaf);
g2=w(nitrogen,carbon,root);

% GRAPHICAL OUTPUT 
figure(1) % Carbon uptake rate & N uptake rate
plot(T,kc,'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(T,kn,'LineWidth',4);
hold on, drawnow
xlabel('Time','FontSize',45)
ylabel('Kc','FontSize',45)
set(gca,'LineWidth',2,'FontSize',45)
legend('K_c 0','K_n 0','K_c 1','K_n 1','K_c 2','K_n 2','K_c 3','K_n 3','K_c 4','K_n 4','K_c 5','K_n 5','Location','Best')
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.65 ;y_width=14.99;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('0c1resp1','-depsc','-loose');

figure(2) % Intermediate C and N
plot(T,carbon,'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(T,nitrogen,'LineWidth',4);
hold on, drawnow
xlabel('Time','FontSize',45)
ylabel('Intermediate concentration','FontSize',45)
set(gca,'LineWidth',2,'FontSize',45)
legend('Carbon 0','Nitrogen 0','Carbon 1','Nitrogen 1','Carbon 2','Nitrogen 2','Carbon 3','Nitrogen 3','Carbon 4','Nitrogen 4','Carbon 5','Nitrogen 5','Location','Best')
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.04 ;y_width=13.12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('0c1resp2','-depsc','-loose');

figure(3) % RGR 
plot(T,rgr1,'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(T,rgr2,'LineWidth',4);
hold on, drawnow
xlabel('Time','FontSize',45)
ylabel('Relative growth rate (RGR)','FontSize',45)
set(gca,'LineWidth',2,'FontSize',45)
legend('Leaf 0','Root 0','Leaf 1','Root 1','Leaf 2','Root 2','Leaf 3','Root 3', 'Leaf 4','Root 4','Leaf 5','Root 5','Location','Best')
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.04 ;y_width=13.12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('0c1resp3','-depsc','-loose');

figure(4) % Leaf and root mass 
plot(T,leaf,'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(T,root,'LineWidth',4);
hold on, drawnow
xlabel('Time','FontSize',45)
ylabel('Mass','FontSize',45)
set(gca,'LineWidth',2,'FontSize',45)
legend('Leaf 0','Root 0','Leaf 1','Root 1','Leaf 2','Root 2','Leaf 3','Root 3', 'Leaf 4','Root 4','Leaf 5','Root 5','Location','Best')
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.04 ;y_width=13.12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('0c1resp4','-depsc','-loose');

figure(5) % Leaf and root mass 
plot(T,SR,'LineWidth',4);
hold on, drawnow 
xlabel('Time','FontSize',45)
ylabel('Shoot:Root','FontSize',45)
set(gca,'LineWidth',2,'FontSize',45)
legend('Leaf 0','Root 0','Leaf 1','Root 1','Leaf 2','Root 2','Leaf 3','Root 3', 'Leaf 4','Root 4','Leaf 5','Root 5','Location','Best')
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.04 ;y_width=13.12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('0c1resp5','-depsc','-loose');

figure(6) % Growth rate  
plot(T,g1,'marker','square','MarkerSize',3,'LineWidth',4);
hold on, drawnow
plot(T,g2,'LineWidth',4);
hold on, drawnow
xlabel('Time','FontSize',45)
ylabel('Growth rate','FontSize',45)
set(gca,'LineWidth',2,'FontSize',45)
legend('Leaf 0','Root 0','Leaf 1','Root 1','Leaf 2','Root 2','Leaf 3','Root 3', 'Leaf 4','Root 4','Leaf 5','Root 5','Location','Best')
 set(gcf, 'PaperUnits', 'inches');
 x_width=15.04 ;y_width=13.12;
 set(gcf, 'PaperPosition', [0 0 x_width y_width]);
print('0c1resp6','-depsc','-loose');
% carbon uptake rate
    function ff=f(x)
        ff=cA.*x./(cx0+x);
    end
%nitrogen uptake rate
    function gg=g(x)
        gg=nA.*x./(nx0+x);
    end
% leaf growth rate
    function hh=h(x,y)
        hh=((v.*x)./(k+x)).*((y.*v)./(y+k));
    end
% root growth rate 
    function dd=d(x,y)
        dd=((v.*x)./(x+k)).*((v.*y)./(y+k));
    end
%leaf growth *mass
    function pp=p(x,y,z)
        pp=((v.*x)./(k+x)).*((y.*v)./(y+k)).*z;
    end
% root mass *mass
        function ww=w(x,y,z)
            ww=((v.*x)./(k+x)).*((y.*v)./(y+k)).*z;
        end
% system of odes
function dydt=thorn(~,y)
    dydt=zeros(4,1);
    dydt(1)=f(y(3))-h(y(1),y(2))*y(3)-d(y(2),y(1)).*y(4)-r1.*(y(3))-r2.*(y(4));
    dydt(2)=g(y(4))-h(y(1),y(2)).*y(3)-d(y(2),y(1)).*y(4);
    dydt(3)=2.*h(y(1),y(2)).*y(3);
    dydt(4)=2.*d(y(2),y(1)).*y(4);
end
end