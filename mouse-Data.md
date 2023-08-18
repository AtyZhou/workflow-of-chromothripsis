# Mouse data bulkRNAseq workflow

This Pipeline is deal with **bulkRNAseq** mouse data start with fastq file to the fusion-gene, SNP and CNV in CC-LY lab

## Fusion-gene
In this Pipeline we use three software: [**fusioncatcher**](https://github.com/ndaniel/fusioncatcher), [**STAR-Fusion**](https://github.com/STAR-Fusion/STAR-Fusion) and [**ericscript**](https://github.com/smsrts/EricScript) to predict fusion-gene

- fusioncatcher
~~~shell
#9022 server 
#create fusioncatcher: conda create -n fusioncatcher fusioncatcher
#go to conda install pathway: cd ~/miniconda3/bin
source activate fusioncatcher
conda deactivate

#fusioncatcher's reference data:
data_source=/mnt/data/user_data/zhaolei/program/yes/envs/fusioncatcher/data/mus_musculus
#output path:
output_path=/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/

#run
fusioncatcher \
-p 20 \
-d $data_source  \
-i /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/1_Sine/ \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/1_Sine/

screen -S fusioncatcher
cd ~/miniconda3/bin
conda activate fusioncatcher
cd /mnt/data/user_data/shuhao/project/RNA eq/chromothripsis/WXD_chromothripsis/result

fusioncatcher \
-p 20 \
-d $data_source  \
-i /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/1_TSP/ \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/1_TSP/

fusioncatcher \
-p 20 \
-d $data_source  \
-i /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/2_Sine/ \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/2_Sine/

fusioncatcher \
-p 20 \
-d $data_source  \
-i /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/2_TSP/ \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/2_TSP/

fusioncatcher \
-p 20 \
-d $data_source  \
-i /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/3_Sine/ \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/3_Sine/

fusioncatcher \
-p 20 \
-d $data_source  \
-i /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/3_TSP/ \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/3_TSP/

~~~
- STAR-Fusion
~~~shell
#9022 server
#install
conda search star-fusion
conda create -n star_fusion star-fusion=1.11.1
conda install -c bioconda star-fusion=1.11.1
conda env list 
cd ~/miniconda3/bin
source activate 
conda activate star_fusion
conda deactivate

cd /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/

vim work 
1_Sine 1_Sine
1_TSP 1_TSP
2_Sine 2_Sine
2_TSP  2_TSP  
3_Sine  3_Sine  
3_TSP 3_TSP

cat work  |while read id;
do 
echo $id
arr=($id)
path=${arr[1]}
sample=${arr[0]}
fq2=${path}'_2.clean.fq.gz'
fq1=${path}'_1.clean.fq.gz'
echo $fq2
echo $fq1
echo $sample
echo $path
STAR-Fusion --CPU 3 \
     --genome_lib_dir /mnt/data/user_data/zlu/download/Mouse_GRCm39_M31_CTAT_lib_Nov092022.plug-n-play/ctat_genome_lib_build_dir \
     --left_fq /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/$path/$fq1 \
     --right_fq /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/$path/$fq2 \
     --output_dir /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/$sample
done

~~~
- ericscript
~~~shell
#9022 server
#install
conda create -n ericscript ericscript
cd ~/miniconda3/bin
conda activate ericscript 
conda deactivate


screen -S ericscript
ericscript.pl \
-p 25 \
-db /mnt/data/user_data/zlu/reference/ericscript/ericscript_db_musmusculus_ensembl84/ \
--refid mus_musculus -name 1_Sine \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/1_Sine \
/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/1_Sine/1_Sine_1.clean.fq.gz /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/1_Sine/1_Sine_2.clean.fq.gz

cat work  |while read id;
do 
echo $id
arr=($id)
path=${arr[1]}
sample=${arr[0]}
fq2=${path}'_2.clean.fq.gz'
fq1=${path}'_1.clean.fq.gz'
echo $fq2
echo $fq1
echo $sample
echo $path
ericscript.pl \
-p 25 \
-db /mnt/data/user_data/zlu/reference/ericscript/ericscript_db_musmusculus_ensembl84/ \
--refid mus_musculus -name $sample \
-o /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/$sample \
/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/$path/$fq1 /mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/00.CleanData/$path/$fq2
done

cat work |while read id;
do 
echo $id
arr=($id)
path=${arr[1]}
echo $path
mv $path ./ericscript
done 
~~~
- analysis
~~~R
library(chimeraviz)
#说明书https://github.com/ndaniel/fusioncatcher/blob/master/doc/manual.md  结果文件以Spanning_pairs（支持融合的读取对计数（还包括多重映射读取））降序排列
Sine_1_fusioncatcher_data<-read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/1_Sine/final-list_candidate-fusion-genes.txt",header = T,sep="\t")
Sine_1_ericscript_data = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/1_Sine/1_Sine.results.filtered.tsv",header = T,sep="\t")
Sine_1_star_fusion_data  = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/1_Sine/star-fusion.fusion_predictions.tsv",header = F,sep="\t")
dim(Sine_1_star_fusion_data)
Sine_1_fu_work<-Sine_1_fusioncatcher_data[,c(1,2)]
mm<-Sine_1_ericscript_data[,c(1,2)]
mm$GeneName1<-as.character(mm$GeneName1)
mm$GeneName2<-as.character(mm$GeneName2)
Sine_1_er_work<-data.frame(nrow=nrow(mm),ncol=2)
for( i in 1:2){
    for(j in 1:nrow(mm)){
    Sine_1_er_work[j,i]<-toupper(mm[j,i])
    }
}
star_pre<-as.character(Sine_1_star_fusion_data[,1])
work<-intersect(Sine_1_fu_work[,1],Sine_1_er_work[,1])
Sine_1_fu_work[which(Sine_1_fu_work[,1]=="SEC61A1"),]
Sine_1_er_work[which(Sine_1_er_work[,1]=="SEC61A1"),]


Sine_1_path_fusioncatcher <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/1_Sine/final-list_candidate-fusion-genes.txt")
Sine_1_path_eriscript <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/1_Sine/1_Sine.results.filtered.tsv")
Sine_1_path_star_fusion <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/1_Sine/star-fusion.fusion_predictions.tsv")

Sine_1<- import_starfusion(Sine_1_path_star_fusion, "mm10")
plot_circle(Sine_1)

library(RCircos)
data(UCSC.Mouse.GRCm38.CytoBandIdeogram)
data(RCircos.Gene.Label.Data)
RCircos.Set.Core.Components(cyto.info=UCSC.Mouse.GRCm38.CytoBandIdeogram,chr.exclude=NULL,tracks.inside=5,tracks.outside=0)
RCircos.Set.Plot.Area()
RCircos.Chromosome.Ideogram.Plot()

refrencegtf<-"/mnt/data/public_data/reference/Mus/Mus_musculus_UCSC/UCSC/mm10/Annotation/Archives/archive-2015-07-17-14-33-26/Genes/genes.gtf"
txdb<-makeTxDbFromGFF(refrencegtf,format="gtf",circ_seqs=character())


Sine_2_fusioncatcher_data<-read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/2_Sine/final-list_candidate-fusion-genes.txt",header = T,sep="\t")
Sine_2_ericscript_data = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/2_Sine/2_Sine.results.filtered.tsv",header = T,sep="\t")
Sine_2_star_fusion_data  = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/2_Sine/star-fusion.fusion_predictions.tsv",header = T,sep="\t")
dim(Sine_2_fusioncatcher_data)
dim(Sine_2_ericscript_data)
dim(Sine_2_star_fusion_data)

Sine_2_fu_work<-Sine_2_fusioncatcher_data[,c(1,2)]
mm<-Sine_2_ericscript_data[,c(1,2)]
mm$GeneName1<-as.character(mm$GeneName1)
mm$GeneName2<-as.character(mm$GeneName2)
Sine_2_er_work<-data.frame(nrow=nrow(mm),ncol=2)
for( i in 1:2){
    for(j in 1:nrow(mm)){
    Sine_2_er_work[j,i]<-toupper(mm[j,i])
    }
}
star_pre<-as.character(Sine_2_star_fusion_data[,1])
work<-intersect(Sine_2_fu_work[,1],Sine_2_er_work[,1])
Sine_2_fu_work[which(Sine_2_fu_work[,1]=="SEC61A1"),]
Sine_2_er_work[which(Sine_2_er_work[,1]=="SEC61A1"),]

Sine_2_path_fusioncatcher <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/2_Sine/final-list_candidate-fusion-genes.txt")
Sine_2_path_eriscript <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/2_Sine/2_Sine.results.filtered.tsv")
Sine_2_path_star_fusion <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/2_Sine/star-fusion.fusion_predictions.tsv")
Sine_2_fusioncatcher <- import_fusioncatcher(Sine_2_path_fusioncatcher, "mm10")
Sine_2_ericscript=import_ericscript(Sine_2_path_eriscript, "mm10")
Sine_2_starfusion<- import_starfusion(Sine_2_path_star_fusion, "mm10")
plot_circle(Sine_2_fusioncatcher)
plot_circle(Sine_2_ericscript)
plot_circle(Sine_2_starfusion)

Sine_3_fusioncatcher_data<-read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/3_Sine/final-list_candidate-fusion-genes.txt",header = T,sep="\t")
Sine_3_ericscript_data = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/3_Sine/3_Sine.results.filtered.tsv",header = T,sep="\t")
Sine_3_star_fusion_data  = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/3_Sine/star-fusion.fusion_predictions.tsv",header = T,sep="\t")
dim(Sine_3_fusioncatcher_data)
dim(Sine_3_ericscript_data)
dim(Sine_3_star_fusion_data)
TSP_1_fusioncatcher_data<-read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/1_TSP/final-list_candidate-fusion-genes.txt",header = T,sep="\t")
TSP_1_ericscript_data = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/1_TSP/1_TSP.results.filtered.tsv",header = T,sep="\t")
TSP_1_star_fusion_data  = read.table("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/1_TSP/star-fusion.fusion_predictions.tsv",header = T,sep="\t")
dim(TSP_1_fusioncatcher_data)
dim(TSP_1_ericscript_data)
dim(TSP_1_star_fusion_data)
TSP_1_path_fusioncatcher <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/fusioncatcher/1_TSP/final-list_candidate-fusion-genes.txt")
TSP_1_path_eriscript <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/ericscript/1_TSP/1_TSP.results.filtered.tsv")
TSP_1_path_star_fusion <- file.path("/mnt/data/user_data/shuhao/project/RNAseq/chromothripsis/WXD_chromothripsis/result/STAR_FUSION/result/1_TSP/star-fusion.fusion_predictions.tsv")
TSP_1_fusioncatcher <- import_fusioncatcher(TSP_1_path_fusioncatcher, "mm10")
TSP_1_ericscript=import_ericscript(TSP_1_path_eriscript, "mm10")
TSP_1_starfusion<- import_starfusion(TSP_1_path_star_fusion, "mm10")
plot_circle(TSP_1_fusioncatcher)
plot_circle(TSP_1_ericscript)
plot_circle(TSP_1_starfusion)
~~~












