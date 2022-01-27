#!/data2/wangb/anaconda2/bin/Rscript
# when dotplot/barplot/..., BP[[1]]@result$Description is sorted by its BP[[1]]@result$Count
rm(list = ls())
options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS=FALSE)


library(DOSE)
library(GO.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(topGO)
library(GSEABase)
library(dplyr)
library(magrittr)
library(clusterProfiler)
library(ggplot2)


enrichIDWithOnt_AndPlot <- function(ids, ont, output_prefix, pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, showCategory=12, OrgDb='org.Hs.eg.db', fontsize=20){
	# Enrichment. 
	result_enrichGO <- enrichGO(gene = ids, OrgDb = OrgDb, ont = ont, 
				pAdjustMethod = pAdjustMethod, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff, readable = TRUE) 
	# result_enrichGO <- as.data.frame(result_enrichGO)
	
	path_output_xls = paste(output_prefix, ont,'xls', sep='.')
	print('output xls to ')
	print(path_output_xls)
	write.table(result_enrichGO@result, path_output_xls, sep='\t')
	
	df = result_enrichGO@result
	print(dim(df))
	df = df[which(df['p.adjust'] < pvalueCutoff),]
	if(dim(df)[1] > 0){
		df = df[which(df['qvalue'] < qvalueCutoff),]
	}
	print(dim(df))

	
	if(nrow(df) == 0){
		print('# Early stop, not plotting')
		return(result_enrichGO)
	}
	
	# print('result_enrichGO')
	# print(result_enrichGO)
	print(paste('ploting', ont, sep=' '))
	
	print('barplot')
	pdf(paste(output_prefix, ont, 'barplot', 'pdf', sep='.'), PIC_WIDTH,PIC_HEIGHT) #, res = NA) 
	print(barplot(result_enrichGO, showCategory=showCategory, font.size=fontsize))
	dev.off()
	
	print('dotplot')
	pdf(paste(output_prefix, ont, 'dotplot', 'pdf', sep='.'), PIC_WIDTH,PIC_HEIGHT) #, res = NA)
	print(dotplot(result_enrichGO, font.size=fontsize))
	dev.off()
	
	if (ont != 'ALL'){
	print('plotGOgraph ')
	pdf(paste(output_prefix, ont, 'plotGOgraph', 'pdf', sep='.'), PIC_WIDTH,PIC_HEIGHT) #, res = NA)
	plotGOgraph(result_enrichGO)
	dev.off()
	}
	# Error.
	# pdf(paste(output_prefix, ont, 'enrichMap', 'pdf', sep='.'), PIC_WIDTH,PIC_HEIGHT)
	# enrichMap(result_enrichGO)
	# dev.off()

	# Error.
	# print('cnetplot')
	# pdf(paste(output_prefix, ont, 'cnetplot', 'pdf', sep='.'), PIC_WIDTH,PIC_HEIGHT)
	# cnetplot(result_enrichGO, categorySize="pvalue", foldChange=ids)
	# dev.off()
	return(result_enrichGO)
	
}

bitrFromDirtySymbolToOtherType <- function(gene_symbols, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db'){
	df_perfect_match <- bitr(gene_symbols, fromType=fromType, toType=toType, OrgDb=OrgDb)
	if (fromType == 'SYMBOL'){
		symbols_not_perfect <- setdiff(gene_symbols, df_perfect_match[, fromType])
		print('symbol cannot map:')
		print(symbols_not_perfect)
		symbols_from_alias = NULL
		tryCatch(
			{
				symbols_from_alias <- bitr(symbols_not_perfect, fromType="ALIAS", toType=fromType, OrgDb=OrgDb)
			}, error = function(e){
				symbols_from_alias <- data.frame(1,2)
				colnames(symbols_from_alias) <- c(fromType, 'ALIAS')
			}
		) 
		print('alias to symbols:')
		print(symbols_from_alias[,fromType]) 
		map_failed <- setdiff(symbols_not_perfect, symbols_from_alias$ALIAS)
		
		print('map failed: ')
		print(map_failed)
		
		df_alias_match <- bitr(symbols_from_alias[,fromType], fromType=fromType, toType=toType, OrgDb=OrgDb)
		print('df_perfect_match')
		print(df_perfect_match)
		print('df_alias_match')
		print(df_alias_match)
		df_to_return <- rbind(df_perfect_match, df_alias_match)
		print('df_converted')
		print(df_to_return)
	}else{
		df_to_return = df_perfect_match
	}
	return(df_to_return)
}

# entriz_list <- bitrFromDirtySymbolToOtherType(gene_symbols)  
kegg_plot <- function(diff_kegg){
	dat<-diff_kegg
	colnames(dat)
	dat$pvalue = -log10(dat$pvalue)
	dat$pvalue=dat$pvalue*dat$group 
	dat=dat[order(dat$pvalue,decreasing =F),]
	g_kegg<- ggplot(dat, aes(x=reorder(Description,order(pvalue, decreasing =T)), y=pvalue, fill=group)) + 
		geom_bar(stat="identity") + 
		scale_fill_gradient(low="blue",high="red",guide = FALSE) + 
		scale_x_discrete(name ="Pathway names") +
		scale_y_continuous(name ="log10P-value") +
		coord_flip() + theme_bw()+
		theme(text = element_text(size=8),plot.title = element_text(hjust = 0.1))+
		ggtitle("Pathway Enrichment") 
}

plotEnrichmentsByGeneSymbols <- function(gene_symbols, prefix, fromType='SYMBOL', toType='ENTREZID', OrgDb='org.Hs.eg.db', do_kegg=T, organism='hsa', showCategory=12, fontsize=20, pvalueCutoff = 0.05, qvalueCutoff = 0.05){
	# gene_ids <- bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
	gene_ids <- bitrFromDirtySymbolToOtherType(gene_symbols, fromType=fromType, toType=toType, OrgDb=OrgDb)
	gene_ids <- gene_ids$ENTREZID
	print('gene ids for GO: ')
	print(gene_ids)
	
	result_groupGO <- groupGO(gene = gene_ids, OrgDb = OrgDb, readable = TRUE)
	pdf(paste(prefix, 'groupGO', 'pdf', sep='.'), PIC_WIDTH, PIC_HEIGHT) #, res = NA)
	print(barplot(result_groupGO, drop=TRUE, showCategory=showCategory))
	dev.off()
	
	enrichment_results <- c()
	
	for(ont in c("BP", "CC", "MF")){
		enrich_result <- enrichIDWithOnt_AndPlot(gene_ids, ont, prefix, OrgDb=OrgDb, fontsize=fontsize, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
		enrichment_results <- c(enrichment_results, enrich_result)
	}
	# return(result_groupGO)
	if (do_kegg){ 
		kk.diff <- enrichKEGG(gene=gene_ids,organism=organism, pvalueCutoff =pvalueCutoff, qvalueCutoff = qvalueCutoff)

		# ggsave(paste0(prefix, '.kegg.dotplot.pdf'),width=10,height=10)#鐐瑰浘
		print(kk.diff)
		barplot(kk.diff)
		ggsave(paste0(prefix, '.kegg.barplot.pdf'),width=10,height=10) 
		
		kegg_diff_dt <- as.data.frame(kk.diff)
		print(kegg_diff_dt)
		if(nrow(kegg_diff_dt) == 0){
			print('# Early stop, not plotting kegg')
			return(enrichment_results)
		}
		
		diff_kegg<-kegg_diff_dt
		diff_kegg$group=-1
		write.table(diff_kegg,paste0(prefix, '.kegg.xls'), sep='\t')
		
		diff_kegg<-kegg_diff_dt[kegg_diff_dt$pvalue<pvalueCutoff,]
		print(diff_kegg)
		if(nrow(diff_kegg) > 0){
			g_kegg<-kegg_plot(diff_kegg) 
			dotplot(kk.diff)
			ggsave(g_kegg,filename = paste0(prefix, '.kegg.pdf'))
		}
		# https://www.jianshu.com/p/79a49da8ea16
		#require(pathview)
		#geneList<-gene_ids#閫夊畾鐨勫熀鍥?
		# pathview(gene.data = "genelist", pathway.id = "hsa04110",
						# species="hsa",limit =list(gene=max(abs(geneList)),cpd=1))
		#pathview(gene.data = "genelist", pathway.id = "hsa04110",
		#				species="hsa",limit=list(gene=100,cpd=1))
	} 
	return(enrichment_results)
}


if(sys.nframe() == 0) {
	args <- commandArgs(T)
	file_in = args[1]

	if (length(args) > 1){
		fromType = args[2]
	}else{
		fromType = 'SYMBOL'
	}

	if (length(args) > 2){
		toType = args[3]
	}else{
		toType = 'ENTREZID'
	} 
	if (length(args) > 3){
		OrgDb = args[4]
	}else{
		OrgDb = "org.Hs.eg.db"
	}
	if (length(args) > 4){
		organism = args[5]
	}else{
		organism = 'hsa'
	}
	if (length(args) > 5){
		pvalueCutoff = as.numeric(args[6])
	}else{ 
		pvalueCutoff = 0.1
	}
	if (length(args) > 6){
		qvalueCutoff = as.numeric(args[7])
	}else{
		qvalueCutoff = 0.1
	}
	
	print('qvalueCutoff:',)

	# file_in should be a text file with one gene symbol each row. 
	file_in = "C://Users//thecat//Desktop//thecat//富集通路分析//基因.txt"

	# PIC_WIDTH = 20
	# PIC_HEIGHT = 20
	showCategory = 12
	# fontsize=30
	scale = 0.5
	# PIC_WIDTH = showCategory * scale + fontsize / 4 
	# PIC_HEIGHT = showCategory * scale
	fontsize=20
	PIC_WIDTH = showCategory * scale + fontsize / 4
	PIC_HEIGHT = showCategory * scale
	
	gene_symbols <- scan(file=file_in, what=" ", sep="\t")
	result_groupGO <- plotEnrichmentsByGeneSymbols(gene_symbols, file_in, fromType=fromType, toType=toType, do_kegg=T, OrgDb=OrgDb, organism=organism, showCategory=showCategory, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff) 
	# result_groupKegg <- plotEnrichmentsByGeneSymbols(gene_symbols, file_in, fromType=fromType, toType=toType, OrgDb=OrgDb, organism=organism, showCategory=showCategory, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff) 
}



