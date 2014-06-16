require("goseq")
require("org.Hs.eg.db")
require("biomaRt")
require("GO.db")
require("pathview")

getGo <- function(de.genes, GO.type) {
  de.genes.fc <- de.genes[,c("gene", "log2.fold_change.")]
  de.genes.fc<-de.genes.fc[!is.infinite(de.genes.fc$log2.fold_change.),]
  de.genes.info <- select(org.Hs.eg.db, keys=as.vector(de.genes.fc$gene),
                        keytype="SYMBOL", columns=c("ENSEMBL", "ENTREZID"))

  # Remove genes without ENSEMBL IDs.
  de.genes.info<-de.genes.info[!is.na(de.genes.info$ENSEMBL),]

  mart<-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  all.genes<-getBM(attributes=c('ensembl_gene_id', 'start_position',
                              'end_position'), mart=mart)
  gene.vector<-as.integer(
                all.genes$ensembl_gene_id%in%de.genes.info$ENSEMBL)
  names(gene.vector) <- all.genes$ensembl_gene_id

  all.genes$length<-all.genes$end_position - all.genes$start_position

  gene.vector.length<-all.genes$length
  names(gene.vector.length)<-all.genes$ensembl_gene_id
  pwf=nullp(gene.vector, bias.data=gene.vector.length)

  GO = goseq(pwf, "hg19", "ensGene", test.cats=GO.type)
  GO$padjust = p.adjust(GO$over_represented_pvalue,
                        method="BH")

  # Get pathway names for significant patways
  GO.sig = GO[GO$padjust<0.05,]
  terms = stack(lapply(mget(GO.sig$category, GOTERM), Term))
  GO.sig$terms = terms$values

  de.genes.fc.path <- merge(de.genes.info, de.genes.fc,
                            by.x="SYMBOL", by.y="gene")
  return(list(GO.sig, de.genes.fc.path))
}
# Write genes in each pathway to separate files
getGoGenes <- function(cat, data, prefix) {
  m = match(xx[[cat]], data$ENTREZID)
  mm = m[!is.na(m)]
  if(length(mm) == 0) { return() }
  d = data.frame(cat, data[mm, ]$ENSEMBL, data[mm,]$ENTREZID,
                 data[mm,]$SYMBOL, data[mm, ]$log2.fold_change.)
  filename = paste(prefix, cat, sep="_")
  write.table(d, filename, sep="\t", row.names=F,
              col.names=F, quote=F)
  return(d)
}
xx = as.list(org.Hs.egGO2ALLEGS)
xx = xx[!is.na(xx)] # remove GO IDs that do not match any gene

data <- read.table('gene_exp.diff', header=T, sep='\t',
                   stringsAsFactor=F)
go.results <- getGo(data[which(data$log2.fold_change. >= 2.0),],
                    "GO:MF")

df = lapply(go.results[[1]]$category,
            getGoGenes, go.results[[2]], "molfunc_up_genes")
write.table(go.results[[1]], 'molfunc_up_terms.txt', sep='\t', quote=F)

go.results <- getGo(data[which(data$log2.fold_change. <= 2.0),],
                    "GO:MF")

df = lapply(go.results[[1]]$category,
            getGoGenes, go.results[[2]], "molfunc_down_genes")
write.table(go.results[[1]], 'molfunc_down_terms.txt', sep='\t', quote=F)

go.results <- getGo(data[which(data$log2.fold_change. >= 2.0),],
                    "GO:BP")
df = lapply(go.results[[1]]$category,
            getGoGenes, go.results[[2]], "biopro_up_genes")
write.table(go.results[[1]], 'biopro_up_terms.txt', sep='\t', quote=F)

go.results <- getGo(data[which(data$log2.fold_change. <= 2.0),],
                    "GO:BP")

df = lapply(go.results[[1]]$category,
            getGoGenes, go.results[[2]], "biopro_down_genes")
write.table(go.results[[1]], 'biopro_down_terms.txt', sep='\t', quote=F)
