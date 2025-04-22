close all;
for i=1:1e3
p1 = plot(x(pI),U(i,pI),'-','linewidth',2); hold on
p2 = plot(x(pI),r(pI),'.','linewidth',2);
%plot(x,qc,'-.','linewidth',2);
ax = gca;
set(ax, 'fontsize', 12);
%legend('$u$', '$r(x)$', '$q_c(x)$', 'interpreter','latex')
xlabel('$x$','interpreter','latex');
axis tight;


I = [find(r>0,1,'first'), find(r>0,1,'last')];
Xs = x(I);
line([Xs(1),Xs(1)], ax.YLim,'linestyle','--','color','r','linewidth',2);
line([Xs(2),Xs(2)], ax.YLim,'linestyle','--','color','r','linewidth',2);

legend([p1,p2],'$u$', '$r(x)$', 'interpreter','latex')

drawnow
pause(0.05)
hold off;
end