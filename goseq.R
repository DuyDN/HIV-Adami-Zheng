require("goseq")
require("org.Hs.eg.db")
require('biomaRt')
require('KEGG.db')
require("pathview")


data <- read.table('gene_exp.diff', header=T, sep='\t',
                   stringsAsFactor=F)
de.genes <- data[which(data$log2.fold_change. >= 2.0 |
                         data$log2.fold_change. <= -2.0),]
de.genes.fc <- de.genes[,c("gene", "log2.fold_change.")]
de.genes.fc<-de.genes.fc[!is.infinite(de.genes.fc$log2.fold_change.),]
de.genes.info <- select(org.Hs.eg.db, keys=as.vector(de.genes.fc$gene),
                      keytype="SYMBOL", columns=c("ENSEMBL", "ENTREZID"))

# Remove genes without ENSEMBL IDs.
de.genes.info<-de.genes.info[!is.na(de.genes.info$ENSEMBL),]

mart<-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
all.genes<-getBM(attributes=c('ensembl_gene_id', 'start_position',
                              'end_position'), mart=mart)
gene.vector<-as.integer(all.genes$ensembl_gene_id%in%de.genes.info$ENSEMBL)
names(gene.vector) <- all.genes$ensembl_gene_id

all.genes$length<-all.genes$end_position - all.genes$start_position

gene.vector.length<-all.genes$length
names(gene.vector.length)<-all.genes$ensembl_gene_id
pwf=nullp(gene.vector, bias.data=gene.vector.length)

KEGG = goseq(pwf, "hg19", "ensGene", test.cats="KEGG")
KEGG$padjust = p.adjust(KEGG$over_represented_pvalue,
                        method="BH")

# Get pathway names for significant patways
KEGG.sig = KEGG[KEGG$padjust<0.05,]

pathway = stack(mget(KEGG.sig$category, KEGGPATHID2NAME))

KEGG.sig$pathway = pathway$values

xx = as.list(org.Hs.egPATH2EG)
xx = xx[!is.na(xx)] # remove KEGG IDs that do not match any gene

de.genes.fc.path <- merge(de.genes.info, de.genes.fc, by.x="SYMBOL",
                          by.y="gene")

# Write genes in each pathway to separate files
# Plot all KEGG pathways using pathview
pathviewPlot <- function(gene.data, cat) {  
  pv.out <- pathview(gene.data=gene.data, pathway.id=cat,
                     species="hsa", out.suffix="pathway", same.layer=F,
                     gene.idtype="ENSEMBL", limit=list(gene=4, cpd=1))
}

get_genes_kegg <- function(cat, data, prefix) {
  m = match(xx[[cat]], data$ENTREZID)
  mm = m[!is.na(m)]
  d = data.frame(cat, data[mm, ]$ENSEMBL, data[mm,]$ENTREZID,
                 data[mm,]$SYMBOL, data[mm, ]$log2.fold_change.)
  gene.data <- d[, 5]
  names(gene.data) <- d[,2]
  pathviewPlot(gene.data, cat)
  filename = paste(prefix, cat, sep="_")
  write.table(d, filename, sep="\t", row.names=F,
              col.names=F, quote=F)
  return(d)
}
get_genes_kegg_no_pathview <- function(cat, data, prefix) {
  m = match(xx[[cat]], data$ENTREZID)
  mm = m[!is.na(m)]
  d = data.frame(cat, data[mm, ]$ENSEMBL, data[mm,]$ENTREZID,
                 data[mm,]$SYMBOL, data[mm, ]$log2.fold_change.)
  filename = paste(prefix, cat, sep="_")
  write.table(d, filename, sep="\t", row.names=F,
              col.names=F, quote=F)
  return(d)
}

df = lapply(KEGG.sig[-4,]$category,
            get_genes_kegg, de.genes.fc.path, "pathway_genes")

df = lapply(KEGG.sig[4,]$category,
            get_genes_kegg_no_pathview,
            de.genes.fc.path, "pathway_genes")
