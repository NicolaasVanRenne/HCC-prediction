#	create .gct file
# more info on .gct format on https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
##############################################################################################################################

# this is to create a .gct file from a gene expression table in R. 
# load your gene expression table, columns are samples, rows are features (genes)

#load gene expression table in the object expr
	expr <- my.gene.expr.table
	
#create gct format
	gct <- matrix("",3+nrow(expr),2+ncol(expr))
	gct[1,1] <- "#1.2"  
	gct[2,1] <- nrow(expr)  
	gct[2,2] <- ncol(expr)  
	gct[3,] <- c("ID","annotation",colnames(expr))
	
	gct[4:(3+nrow(expr)), 3:(2+ncol(expr))] <- as.matrix(apply(expr, 2, as.character))
	gct[4:(3+nrow(expr)), 1] <- rownames(expr)
	gct[4:(3+nrow(expr)), 2] <- "na"  
	
#save file
	write.table(gct, file="my_gene_expr_table.gct", quote=F, sep="\t", row.names=F, col.names=F) 
				
