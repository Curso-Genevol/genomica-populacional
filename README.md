Pratica Aula 1 Genômica Populacional
================
Kelly Nunes
02/10/2020

## Preparação para a Prática

# Antes de inciar a prática, vamos carregar dois pacotes de R

``` r
library("gdsfmt")
library("SNPRelate")
```

## Início da Prática

Primeiro vamos indicar em qual pasta os arquivos com os datasets formam
salvos:

Na janela direita abaixo, clique na aba **“Files”**, procure pelo
diretório onde os arquivos para essa aula foram salvos. Em seguida,
clique no icone **“More”** e clique na opção **“Set As Working
Directory”**.

# Explorando os arquivos

Para essa prática usaremos dados públicos, que contem genótipos para
SNPs autossômicos em amostras Nativo Americanas (NAM) do Painel de
Diversidade do Genoma Humano (HGDP) e em quatro populações do HapMap:
Yoruba - Africa (YRI); Utah - Europa (CEU), Mexicanos residentes em Los
Angeles (MXL) e Afro-Americanos - EUA (ASW)

  - Qual o número de indivíduos nesse estudo?

<!-- end list -->

``` r
IND<-read.table(file="Samples.fam",sep=" ", header=FALSE,na="NA")
head(IND)
dim(IND)
unique(IND$V1)
```

  - Qual o número de SNPs

<!-- end list -->

``` r
SNPs<-read.table(file="Samples.bim",sep="\t", header=FALSE,na="NA")
head(SNPs)
dim(SNPs)
```

  - Qual o número de amostras em cada população?

<!-- end list -->

``` r
POPINFO=read.table(file="Population_Sample_Info.txt",header=TRUE)
table(POPINFO$Population)
```

# Análise para verificar estruturação populacional

Vamos iniciar a análise de estrutura populacional explorando com PCA.

Primeiro é necessário converter os arquivos para o formato gds. O
formato gds é um formato criado para o uso do pacote de R SNPRELATE.
Esse formato visa organizar os dados de modo a otimizar o tempo de
excecução das análises.

``` r
bedfile<-"Samples.bed" 
famfile<-"Samples.fam" 
bimfile<-"Samples.bim"

snpgdsBED2GDS(bedfile, famfile, bimfile, "Geno.gds")
```

Em seguida, podemos verificar como é o arquivo no formato gds

``` r
genofile <- snpgdsOpen("Geno.gds")
head(genofile)

head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))
```

Por fim, vamos realizar a análise de PCA e discutir o resultado através
de um gráfico

``` r
pca <- snpgdsPCA(genofile)

sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
population=POPINFO$Population

tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(population)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)


plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1", main="PCA using all SNPs")
legend("topleft", legend=levels(tab$pop), pch="o", col=1:(nlevels(tab$pop)))
```

  - Quanto cada PC explica da variância genética observada?

<!-- end list -->

``` r
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))
lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
pairs(pca$eigenvect[,1:4], col=tab$pop, labels=lbls)
```

Em nossa amostra, temos duas populações do continente Americanos que
foram formadas pelo processo de miscigenação: Afro-Americanas e
Mexicanas.

Quais os componentes de ancestralidade de cada uma dessas populações?

  - Vamos, com base nas análises de PCA estimar a proporção de
    componente de ancestralidade Nativo Americano e Europeu nos
    indivíduos Mexicanos.

<!-- end list -->

``` r
avgCEU2=mean(pca2$eigenvect[population=="CEU",2])
avgNAM2=mean(pca2$eigenvect[population=="NAM",2])

MXLadmix=(pca2$eigenvect[population=="MXL",2]-avgNAM2)/(avgCEU2-avgNAM2)
```

  - Por fim, vamos visualizar a ancestralidade estimada em gráfico de
    barras

<!-- end list -->

``` r
tab2=cbind(MXLadmix,1-MXLadmix)
myorder=order(MXLadmix)
temp=t(as.matrix(tab2[myorder,]))


barplot(temp, col=c("blue","green"),xlab="Individual ", ylab="Ancestry", border=NA,axisnames=FALSE,main="Ancestry of MXL",ylim=c(0,1))
legend("bottomright", c("European","Native American"), lwd=4, col=c("blue","green"), bg="white",cex=0.85)
```

O que podemos concluir a partir desses gráficos?
