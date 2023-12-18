# CIRI-deep
- CIRI-deep is a deep-learning model used to predict differentially spliced circRNAs between two biological samples using totalRNA sequencing data. 
- An adapted version of CIRI-deep, CIRI-deepA, was trained for poly(A) selected RNA-seq data.

# Installation
The CIRI-deep model was constructed based on Keras. The `environment.yaml` was provided and the dependencies can be installed as the follow:
```
git clone https://github.com/gyjames/CIRIdeep.git
cd CIRIdeep
conda env create -n CIRIdeep -f ./environment.yaml
conda activate CIRIdeep
```

# Usage
The main program `CIRIdeep.py` can be used to predict differentially spliced circRNAs with CIRIdeep or CIRIdeep(A) or train your own model.

## Predict

**Prediction with CIRIdeep using total RNA-seq data**

CIRIdeep provides probability of given circRNAs being differentially spliced between any of two samples. When predict with CIRIdeep, expression value of 1499 RBPs (listed in `./demo/RBPmax_totalRNA.tsv`) and splicing amount (derived from SAM alignment files) in both samples are needed. The order of RBP expression of each sample should keep exactly the same with `RBP max value file`. We recommend to process raw total RNA-seq fastq files with `CIRIquant`, which provides junction ratio of each circRNA and expression value of each gene in a one-stop manual. SAM files generated with BWA is recommended when producing splicing amount values.

```
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax_totalRNA.tsv -seqFeature ./demo/cisfeature.tsv -splicing_max ./demo/splicingamount_max.tsv -predict_list ./demo/predict_list.txt -model_path ./models/CIRIdeep.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_total -splicing_dir ./demo/splicingamount
```

Several files are needed for prediction.

`-geneExp_absmax` This file contains maximum value of 1499 RBP expression value (TPM) across the training datasets used for normalization. 

`-seqFeature` This file contains normalized cis features of circRNAs to be predicted. A table containing cis features of 71459 circRNAs has been constructed.

`-splicing_max` This file contains maximum value of splicing amount of each circRNA across the training datasets used for normalization.

`-predict_list` This file is comprised of two columns. The first column contains the name of sample pairs seperated by `_`. The second column contains the path to files containing circRNA to be predicted.
CircRNAs are given as coodination on `hg19` genome, like `chr10:102683732|102685776`.

`-model_path` We have provided fully trained CIRIdeep model for using.

`-outdir` Directory to output prediction result.

`-RBP_dir` Directory containing the RBP expression value in TPM of samples to be predicted.

`-splicing_dir` Directory containing the splicing amount of circRNAs to be predicted in each sample. We have provided a basic script `script_splicingamount.py` to produce splicing amount in samples.

**Prediction with CIRIdeep(A) using poly(A) selected RNA-seq data**

CIRIdeep(A) gives three probabilities indicating the circRNA being unchanged, having higher junction ratio in sample A or having higher junction ratio in sample B, which sum to one.
Order of samples (A, B) is the same with sample pair name given in  `predict list file`.
As in some cases, like in scRNA-seq or spatial transcriptomics data, only gene expression matrix is provided, splicing amount is not needed in CIRIdeep(A) any more.

```
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax_polyA.tsv -seqFeature ./demo/cisfeature.tsv -predict_list ./demo/predict_list.txt -model_path ./models/CIRIdeepA.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_polyA --CIRIdeepA
```
`--CIRIdeepA` When predict using CIRIdeepA, this parameter is needed.

Basically, the input files are similar to CIRIdeep, excluding splicing amount related files. **Notably**, the `RBP max value file` file is different from that used in CIRIdeep and all the expression values should be derived from poly(A) selected RNA-seq data. Still, when using CIRIdeep(A), the order of RBP expression of each sample should keep exactly the same with `RBP max value file`.

**Generation of input files**

Here we gave necessary instructions for generating the input files from different datasets.

**RBP expression of total RNA-seq data**
There are two columns in RBP expression level file, the first column identify gene symbols and the second column gives expression level of the RBP in TPM. The order of genes should keep exactly the same with `demo/RBPmax_totalRNA.tsv`.

| Gene Name | TPM |
|------|-------|
|A1CF|12.5|
|AAR2|23.9|

**RBP expressin of poly(A) RNA-seq data**

The format is as same as the RBP expression file used in total RNA-seq data. The order of genes should keep exactly the same with `demo/RBPmax_polyA.tsv`.

**RBP expression of single-cell RNA-seq data**

When analyzing differentially spliced circRNA between cell clusters, the mean value of RBP expression level in CPM or TPM was used. The order of genes should keep exactly the same with `demo/RBPmax_polyA.tsv`

**RBP expression of spatial transcriptome data**

We recommend to perform imputation step before extracting expression level of RBPs. Tangram, gimVI and SpaGE were greate choices. After imputation, the gene expression value should be normalized as: $$Exp^i = Exp_{imputed}^i / \Sigma Exp_{imputed}^i*scalefactor$$

We used 300,000 as scale factor here. The order of genes should keep exactly the same with `demo/RBPmax_polyA.tsv`



## Train

**CIRIdeep training**

```
python CIRIdeep.py train -geneExp_absmax /path/to/file -seqFeature /path/to/file -splicing_max /path/to/file -outdir /out/path -RBP_dir /RBP/path -splicing_dir /splicing/path
```
Hyperparameters are given in `config.py`. `config.py` must be under the same directory with `CIRIdeep.py`. Resources are waiting to be loaded...

**CIRIdeep(A) training**

```
python CIRIdeep.py train -geneExp_absmax /path/to/file -seqFeature /path/to/file -outdir /out/path -RBP_dir /RBP/path --CIRIdeepA
```

## Contact
Zihan Zhou. zhouzihan2018m@big.ac.cn

Please open an issue if you find bugs.

