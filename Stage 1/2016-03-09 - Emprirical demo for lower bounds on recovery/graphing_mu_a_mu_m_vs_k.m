k = 5:100;
val = 1-(k/2).*(1-sqrt((k-2)./(k-1)));
plot(k,val,'r-o','Linewidth',2);
set(gca,'fontsize',20)
xlabel('Sparsity, \it{k}','FontSize',28,'fontweight','bold');
ylabel('\mu_a / \mu_m','FontSize',28,'fontweight','bold');