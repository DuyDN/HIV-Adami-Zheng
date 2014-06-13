require("goseq")
require("org.Hs.eg.db")
require('biomaRt')
require('KEGG.db')

data <- read.table('gene_exp.diff', header=T, sep='\t',
                   stringsAsFactor=F)
de.genes <- data[which(data$significant == 'yes'),]
de.genes.fc <- de.genes[,c("gene", "log2.fold_change.")]
rownames(de.genes.fc) <- de.genes$gene
gene.info <- select(org.Hs.eg.db, keys=as.vector(data.sig$gene),
                      keytype="SYMBOL", columns=c("ENSEMBL", "ENTREZID"))
all.genes <- select(org.Hs.eg.db, keys=keys(org.Hs.eg.db),
                       columns="ENSEMBL")
gene.vector<-as.integer(all.genes$ENSEMBL%in%gene.info$ENSEMBL)
names(gene.vector) <- gene.info$ENSEMBL
mart<-useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

mart<-useMart(biomart="ensembl", dataset="ggallus_gene_ensembl")
all.genes<-getBM(attributes=c('ensembl_gene_id', 'start_position',
                             'end_position'), mart=mart)

all.genes$length<-all.genes$end_position - all.genes$start_position

gene.vector<-as.integer(all.genes%in%gene.info$ENSEMBL)
names(gene.vector)<-all.genes

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

de.genes.fc.path <- merge(gene.info, de.genes.fc, by.x="SYMBOL",
                          by.y="gene")

# Write genes in each pathway to separate files
get_genes_kegg = function(cat, data, prefix)
{
  m = match(xx[[cat]], data$ENTREZID)
  mm = m[!is.na(m)]
  d = data.frame(cat, data[mm, ]$ENSEMBL, data[mm,]$ENTREZID,
                 data[mm,]$SYMBOL, data[mm, ]$log2.fold_change.)
                 
  filename = paste(prefix, cat, sep="_")
  write.table(d, filename, sep="\t", row.names=F,
              col.names=F, quote=F)
  return(d)
}

df = lapply(KEGG.sig$category,
            get_genes_kegg, de.genes.fc.path,
            "pathway_genes")
