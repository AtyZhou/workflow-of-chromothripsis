# Mouse data workflow

This Pipeline is deal with mouse data start from fastq file to the fusion-gene SNP,and CNV in CC-LY lab

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







