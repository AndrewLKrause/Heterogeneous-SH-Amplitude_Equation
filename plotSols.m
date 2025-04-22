close all
h = figure;
h.Position(1) = 0;
h.Position(2) = 0;
h.Position(3) = 1.5*h.Position(3); h.Position(4)=1.5*h.Position(4);
%pI = find(x< 0.51 & x > 0.49);
pI = x>-1;
p1 = plot(x(pI),U(end,pI),'-','linewidth',1); hold on
%p2 = plot(x(pI),r(pI),'.','linewidth',2);
%plot(x,qc,'-.','linewidth',2);
ax = gca;
set(ax, 'fontsize', 16);
%legend('$u$', '$r(x)$', '$q_c(x)$', 'interpreter','latex')
xlabel('$x$','interpreter','latex');

axis tight manual;

Xs = AllZeros(rf, 0, 1, N);

for i=1:length(Xs)
    line([Xs(i),Xs(i)], ax.YLim,'linestyle','--','color','r','linewidth',1);
end
if(subcritical)
    for i=1:length(Xs2)
        line([Xs2(i),Xs2(i)], ax.YLim,'linestyle','--','color','g','linewidth',1);
    end
    for i=1:length(Xs3)
        line([Xs3(i),Xs3(i)], ax.YLim,'linestyle','--','color','b','linewidth',1);
    end
end
%legend([p1,p2],'$u$', '$r(x)$', 'interpreter','latex')
%legend([p1], '$u$', 'interpreter','latex')
