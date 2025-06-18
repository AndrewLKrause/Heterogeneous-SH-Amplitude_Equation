
close all;
f1=figure;

maxU = max(max(U)); minU = min(min(U));
for i=1:length(T)
    if ~ishghandle(f1)
        break
    end
    hold off
    plot(x,U(i,:),'linewidth',2); hold on
    plot(x, A(i,:),'--k','linewidth',2)
    plot(x, -A(i,:),'--k','linewidth',2)
    set(gca,'fontsize',24);
    axis tight;
    set(gca,'YLim',[minU,maxU]);
    title(['$t = ',num2str(T(i)),'$'],'interpreter','latex')
    xlabel('$x$','interpreter','latex');
    pause(2/length(T));
end

