###############################
### Prática Curso de Férias ###
###############################
   
getwd()
setwd("/Users/Kelly/Documents/Class/_Curso_Ferias/Pratica/")

## Preparação para a Prática 

#Antes de inciar a prática, vamos carregar dois pacotes de R
library("gdsfmt")
library("SNPRelate")


## Início da Prática

# Primeiro vamos indicar em qual pasta os arquivos com os datasets formam salvos:
        
# Na janela direita abaixo, clique na aba **"Files"**, procure pelo diretório onde os arquivos para essa aula foram salvos. Em seguida, clique no icone **"More"** e clique na opção **"Set As Working Directory"**.


## Explorando os arquivos ##

        
# Para essa prática usaremos dados públicos, que contem genótipos para SNPs autossômicos em amostras Nativo Americanas (NAM) do Painel de Diversidade do Genoma Humano (HGDP) e em quatro populações do HapMap: Yoruba - Africa (YRI); Utah - Europa (CEU), Mexicanos residentes em Los Angeles (MXL) e Afro-Americanos - EUA (ASW) 


# - Qual o número de indivíduos nesse estudo?

IND<-read.table(file="Samples.fam",sep=" ", header=FALSE,na="NA")
head(IND)
dim(IND)


# - Qual o número de SNPs

SNPs<-read.table(file="Samples.bim",sep="\t", header=FALSE,na="NA")
head(SNPs)
dim(SNPs)


# - Qual o número de amostras em cada população?

POPINFO=read.table(file="Population_Sample_Info.txt",header=TRUE)
table(POPINFO$Population)


# Como é possível identificar estruturação populacional? ##

        
# Vamos iniciar a análise de estrutura populacional explorando com PCA (Análise de componente Principal - Principal Component Analyses).


# A abordagem de PCA é amplamente utilizada no campo da genômica populacional para inferir estruturação populacional a partir de dados genômicos em larga escala. O PCA irá sumarizar a variação genômica observada entre os individuos. Quanto maior o número de diferenças genéticas entre indivíduos, mais distante eles estarão no PCA. Por outro lado quanto maior for o compartilhamento de variantes genéticas mais próximos os indivíduos estarão.


# Para realizar a análise de PCA, primeiro é necessário converter os arquivos para o formato gds. O formato gds é um formato criado para o uso do pacote de R SNPRELATE. Esse formato visa organizar os dados genômicos de modo a otimizar o tempo de excecução das análises.

bedfile<-"Samples.bed" 
famfile<-"Samples.fam" 
bimfile<-"Samples.bim"

snpgdsBED2GDS(bedfile, famfile, bimfile, "Geno.gds")



# Em seguida, podemos verificar como é o arquivo no formato gds

genofile <- snpgdsOpen("Geno.gds")
head(genofile)

head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))


# Por fim, vamos realizar a análise de PCA e discutir o resultado através de um gráfico

pca <- snpgdsPCA(genofile)

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
population=POPINFO$Population

tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(population)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)


plot(tab$EV2, tab$EV1, col=rainbow(5)[as.integer(tab$pop)], xlab="PC2", ylab="PC1", main="PCA")
legend("topleft", legend=levels(tab$pop), pch="o", col=rainbow(5)[1:(nlevels(tab$pop))], cex=0.7, pt.cex = 0.7)


# - Quanto cada PC explica da variância genética observada?
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=rainbow(5)[tab$pop], labels=lbls)


# - É possível usar apenas um subconjuntosd dos SNPs para estimar os PCs?


# Vamos fazer um filtro excluíndo SNPs em desequilíbrio de ligação.
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

snpset.id <- unlist(snpset)

# Agora, vamos estimar o PCA usando apenas esse subconjunto de SNPS.

pca2 <- snpgdsPCA(genofile, snp.id=snpset.id)

tab2 <- data.frame(sample.id = pca2$sample.id,
                   pop = factor(population)[match(pca2$sample.id, sample.id)],
                   EV1 = pca2$eigenvect[,1],    # the first eigenvector
                   EV2 = pca2$eigenvect[,2],    # the second eigenvector
                   stringsAsFactors = FALSE)

head(tab2)

plot(tab2$EV2, tab2$EV1, col=rainbow(5)[as.integer(tab2$pop)], xlab="PC2", ylab="PC1", main="PCA")
legend("bottomleft", legend=levels(tab2$pop), pch="o", col=rainbow(5)[1:(nlevels(tab2$pop))], cex=0.7, pt.cex = 0.7)



## Em nossa amostra, temos duas populações do continente Americanos que foram formadas pelo processo de miscigenação: Afro-Americanas e Mexicanas.

# Quais os componentes de ancestralidade de cada uma dessas populações?
   
# Vamos, com base nas análises de PCA estimar a proporção de componente de ancestralidade Africano e Europeu nos indivíduos Afro-Americanos (ASW).

avgYRI <-mean(pca2$eigenvect[population=="YRI",2])
avgCEU <-mean(pca2$eigenvect[population=="CEU",2])
ASWadmix<-(pca2$eigenvect[population=="ASW",2]-avgCEU)/(avgYRI-avgCEU)

# Ancestralidade média 
AFR_anc<-mean(ASWadmix)
AFR_anc
EUR_anc<-1-AFR_anc
EUR_anc

# - Por fim, vamos visualizar a ancestralidade estimada em gráfico de barras
tab<-cbind(ASWadmix,1-ASWadmix)
myorder<-order(ASWadmix)
temp<-t(as.matrix(tab[myorder,]))

# Visualização Grafica
barplot(temp, col=c("purple","yellow"), xlab="Individual ", ylab="Ancestry", border=NA, axisnames=FALSE, main="Ancestry of MXL",ylim=c(0,1), xpd = FALSE)
legend("bottomright", c("African","European"), lwd=4, col=c("purple","yellow"), bg="white",cex=0.85)



# Vamos, com base nas análises de PCA estimar a proporção de componente de ancestralidade Nativo Americano e Europeu nos indivíduos Mexicanos.
avgCEU=mean(pca2$eigenvect[population=="CEU",2])
avgNAM2=mean(pca2$eigenvect[population=="NAM",2])

MXLadmix=(pca2$eigenvect[population=="MXL",2]-avgNAM2)/(avgCEU-avgNAM2)


# - Por fim, vamos visualizar a ancestralidade estimada em gráfico de barras
tab2=cbind(MXLadmix,1-MXLadmix)
myorder=order(MXLadmix)
temp=t(as.matrix(tab2[myorder,]))

# Ancestralidade média
EUR_anc2<-mean(MXLadmix)
EUR_anc2
NAM_anc<-1-EUR_anc2
NAM_anc

# Visualização gráfica
barplot(temp, col=c("yellow","blue"), xlab="Individual ", ylab="Ancestry", border=NA, axisnames=FALSE, main="Ancestry of MXL",ylim=c(0,1), xpd = FALSE)
legend("bottomright", c("European","Native American"), lwd=4, col=c("yellow","blue"), bg="white",cex=0.85)


# O que podemos concluir a partir desses gráficos?