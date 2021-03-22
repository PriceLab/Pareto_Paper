function signif_genes = Show_Significant_Genes(data,arcOrig,n,gene_list)
    arcnum=size(arcOrig,1); 
    signif_genes=[];
    for i=1:arcnum
        ra(i,:)=abs(arcOrig(i,:)-mean(data));
        [~,idx] = sort(ra(i,:),'descend');
        m = mean(data);
        signif_genes = [signif_genes, gene_list(idx(1:n)), num2cell(arcOrig(i,idx(1:n))-m(idx(1:n)))'];
        titl{i} = ['arc ',num2str(i)];
    end
