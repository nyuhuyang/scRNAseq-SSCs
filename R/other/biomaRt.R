library(biomaRt)
listMarts()
ensembl=useMart("ensembl")
listDatasets(ensembl)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
filters = listFilters(ensembl)
filters[1:20,]
attributes = listAttributes(ensembl)
attributes[grepl("Gene name",attributes$description),]
attributes[grepl("UniProtKB",attributes$description),]
selected_attri <- c('uniprot_gn',
                    'external_gene_name')
lapply(selected_attri,function(x) attributes[grepl(x,attributes$name),])
F_Box_Protein_IDs = read.csv("./doc/F_Box_Protein_IDs.csv",header = T)
biomart_results <- getBM(attributes=selected_attri, 
                         filters= "uniprot_gn",
                         values=F_Box_Protein_IDs$gene_id,
                         mart = ensembl)
colnames(biomart_results)[colnames(biomart_results) == "external_gene_name"] = "Symbol"
write.csv(biomart_results,"./output/Paula20180813/F_Box_Protein_IDs.csv")
