function h = SortAndBarPlot(signif_genes)

for i=1:size(signif_genes,2)/2
    [~,idx] = sort(cell2mat(signif_genes(:,i*2)),'descend');  
    h(i) = figure('PaperOrientation','landscape','Units','characters','Position',[50 50 180 20]);
    bar(cell2mat(signif_genes(idx,i*2)))
    set(gca,'xticklabel',strcat(signif_genes(idx,2*i-1),{' '}),'xtick',1:size(signif_genes,1),'xlim',[0,size(signif_genes,1)+1],...
        'XGrid','on','FontSize',8)
    %set(gca,'xticklabel',gene_list,'xtick',1:length(gene_list),'XGrid','on','FontSize',8)
    %rotateXLabels( gca(), 90 )
    ax = gca;
    ax.XTickLabelRotation = 90;
    title(['arc ', num2str(i)])
end


