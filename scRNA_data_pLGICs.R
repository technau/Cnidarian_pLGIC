# pLGIC paper
library(easypackages)
libraries("Matrix", "readxl","RColorBrewer",
          'patchwork','dplyr','viridisLite','ggplot2','pals','Seurat')

# run this section for set-up
  ## load genes 
  load(file='Genes.Nv2.RData') #gene annotations and colour palettes
  load(file='Alldata.Nv2.publish.Robj') #published expression atlas from sea-anemone-atlas.cells.ucsc.edu, Nv2 Atlas

#* Nv2 models corresponding to gene set:
GLR.Nv2=c('NV2.3220','NV2.24331','NV2.3533','NV2.9790','NV2.9791','NV2.8590','NV2.18184','NV2.19891','NV2.22573','NV2.18182','NV2.2026','NV2.3907','NV2.22458')
GABA.a.Nv2.full = c('NV2.10304','NV2.10497','NV2.10940','NV2.11052','NV2.11913','NV2.11917','NV2.1265','NV2.12693','NV2.13369','NV2.13535','NV2.14387','NV2.18183','NV2.18185','NV2.18186','NV2.18574','NV2.19245','NV2.23093','NV2.235','NV2.238','NV2.25664','NV2.371','NV2.7146','NV2.7411','NV2.9691','NV2.9708','NV2.9943','NV2.8594','NV2.22460')
iGluR <- c("NV2.10569", "NV2.11906", "NV2.12584", "NV2.12681", "NV2.13953", "NV2.14488", "NV2.14489", "NV2.14581", "NV2.14687", "NV2.16901","NV2.16902", "NV2.16903", "NV2.17099", "NV2.25060", "NV2.3976", "NV2.4003", "NV2.4129", "NV2.4265", "NV2.5477", "NV2.6566","NV2.6567", "NV2.6572", "NV2.7310", "NV2.8059")

# update names for the current analysis:
new.name = c('Pou4','pLGIC-1:GluCl-1','pLGIC-2:GluCl-2','pLGIC-3','pLGIC-4','pLGIC-5','pLGIC-6:GluCl-3','pLGIC-7')
gene_short_name = c("BRN3-like-1","GBRB2-like-11", "GLRB-like-1", "GBRB4-like-1", "GBRB-like-1", "GLRA3-like-2", "GBRG3-like-1", "GBRA5-like-1")

# subset out the cloned members for future plotting
Gaba.Glu.cloned <- as.data.frame(cbind(new.name, gene_short_name))
genes$gene_short_name[match(Gaba.Glu.cloned$gene_short_name,genes$gene_short_name)]=Gaba.Glu.cloned$new.name
genes$name=paste(genes$geneID,genes$gene_short_name,sep=' : ')

# generate the three lists according to the dataset names:
GABAA=sort(genes$name[match(GABA.a.Nv2.full,genes$geneID)])
GLR= genes$name[match(GLR.Nv2,genes$geneID)]
iGluR = genes$name[match(goi_iGluR,genes$geneID)]

cloned=genes$name[match(Gaba.Glu.cloned$new.name[2:8],genes$gene_short_name)]
cloned = cloned[c(1,2,6,3,4,5,7)]

goi=unique(c(cloned,GLR,GABAA,iGluR))

#define temporary data to preserve the published dataset:
temp=Alldata.Nv2

#update the object names
rownames(temp@assays$RNA@counts)=genes$name
rownames(temp@assays$RNA@data)=genes$name
rownames(temp@assays$RNA@meta.features)=genes$name
temp@assays$RNA@var.features=genes$name

# scale all the goi's by library:
temp<-ScaleData(temp,features = c(goi,'NV2.14459 : FOXL2','NV2.2966 : Pou4','NV2.252 : Elav1', 'NV2.5271 : MELC4'),split.by = 'orig.ident')
modules=list(iGluR,c(GLR,GABAA))
names(modules)=c('iGluR','pLGIC')
temp<-AddModuleScore(temp,modules)

# generate smaller object for figure generation
data1=subset(temp,features = c(goi,'NV2.14459 : FOXL2','NV2.2966 : Pou4','NV2.252 : Elav1', 'NV2.5271 : MELC4'))
data1=DietSeurat(data1,scale.data = T,dimreducs = 'umap')

# merge some Identity clusters for visualization:
data1=SetIdent(data1,value = 'IDs')
data1$reduced=data1@active.ident
levels(data1$reduced)=c("pSC|PGCs", "pSC|PGCs", "neurogladular", "gland.mucous", "unchar.immune", "neurogladular", "nematocyte", "mat.nematocyte", "endoderm", 
                        "ectoderm", "ectoderm", "retractor muscle", "mesoderm", 
                        "mesoderm")
data1$reduced=droplevels(data1$reduced)

# select full identity clusters to keep for display:
data1<-SetIdent(data1,value = 'ID.separate')
coiTR=WhichCells(data1,ident=levels(data1)[c(73,74)]) #merge Tentacle retractor
coiMR=WhichCells(data1,ident=levels(data1)[c(72,75)]) #merge mesentary retractor
coiSprioM=WhichCells(data1,ident=levels(data1)[c(69)]) #keep mature spirocytes separate
coiSprio=WhichCells(data1,ident=levels(data1)[c(68)]) # keep spirocytes separate
data1$ID.separate=as.character(data1$ID.separate)
data1$ID.separate[coiTR]='RM.TR'
data1$ID.separate[coiMR]='RM.MR'
data1$ID.separate[coiSprioM]='mat.spirocyte'
data1$ID.separate[coiSprio]='spirocyte'

data1$ID.separate = as.factor(data1$ID.separate)
data1<-SetIdent(data1,value = 'ID.separate')

data1<-SetIdent(data1,value = 'reduced')

# select the clusters of interest for plotting:
coi=WhichCells(SetIdent(data1,value = 'ID.separate'),idents = levels(SetIdent(data1,value = 'ID.separate'))[c(45,52:57,61:72,74:81,84:90,102:104)])
data1$ID.separate.reduced=data1$reduced
data1$ID.separate.reduced=as.character(data1$ID.separate.reduced)
data1$ID.separate.reduced[coi]=as.character(data1$ID.separate[coi])
data1$ID.separate.reduced=as.factor(data1$ID.separate.reduced)
data1<-SetIdent(data1,value = 'ID.separate.reduced')
x=setdiff(c(1:length(levels(data1))),c(7,4,45,5,9:34,43:44,35:41))
data1@active.ident=factor(data1@active.ident,levels(data1)[c(7,4,45,5,9:34,43:44,35:41,x)])

data1@active.ident<-droplevels(data1@active.ident)
levels(data1)

###Fig7a ----
#subset to the clusters of interest:
data2<-subset(data1,idents = levels(data1)[c(1:32)])
levels(data2)
DotPlot(data2,features=(c(cloned[c(1,2,3)],'NV2.252 : Elav1','NV2.2966 : Pou4', 'NV2.5271 : MELC4')),cols=c('grey','darkred'),col.min = 0, scale.by = 'size',
        #scale.max = 100,scale.min = 0,
        dot.min = 0.01
)&
  RotatedAxis()&FontSize(10,12) &theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&coord_flip()

### Figure S2 ----
#plot all putative Glu receptors across all published single cell data identities:

goi.exp = goi[which(rowSums(temp@assays$RNA@counts[goi,])>=3)] #these all have at least 3 reads...
DotPlot(temp,features=(goi.exp),group.by = 'ID.separate',cols=c('grey','darkred'),col.min = 0, dot.min = 0.01,scale.by = 'size')&
  RotatedAxis()&FontSize(8,10) &ggplot2::theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&labs(title='Expression of all putative GABA-receptors', subtitle ='>3 reads detected across full dataset')

### Figure Sx ----
#filter out non expressed iGlus:
x=DotPlot(data2,features=iGluR)
filtered=x$data$features.plot[x$data$pct.exp>=10]
iGluR.f=unique(as.character(filtered))
#plot all cloned:
p1=DotPlot(data2,features=(c(cloned)),cols=c('grey','darkred'),col.min = 0, scale.by = 'size',
           scale.max = 100,scale.min = 0,
)&
  RotatedAxis() &theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&coord_flip() &theme(axis.text.x= element_blank(),axis.title.x = element_blank())
#Plot filtered iGluRs:
p2=DotPlot(data2,features=(c(iGluR.f)),cols=c('grey','slateblue'),col.min = 0, scale.by = 'size',
        scale.max = 100,scale.min = 0,
)+FontSize(10,12)&
  RotatedAxis() &theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&coord_flip()&plot_layout(ncol=1)

p1/p2

#plot the modules scores for both types: iglu vs plgic:
DotPlot(data2,features=c('Cluster1','Cluster2'),cols=c('grey','slateblue'),col.min = 0, scale.by = 'size',
           scale.max = 100,scale.min = 0,
)+FontSize(10,12)&
  RotatedAxis() &theme(panel.grid = element_line(size=0.002,colour = 'grey90'),legend.title = element_text(size = 8),legend.text = element_text(size=6))&coord_flip()&plot_layout(ncol=1)
#not used.
