# CIRI-deep
- CIRI-deep is a deep-learning model used to predict differentially spliced circRNAs between two biological samples using totalRNA sequencing data. 
- An adapted version of CIRI-deep, CIRI-deepA, was trained for poly(A) selected RNA-seq data.

# Installation
The CIRI-deep model was constructed based on Keras. The `requirements.txt` was provided and the dependencies can be installed as the follow:
```
git clone https://github.com/gyjames/CIRIdeep.git
cd CIRIdeep
conda env create --name CIRIdeep --file ./requirements.txt
```

# Usage
The main program `CIRIdeep.py` can be used to predict differentially spliced circRNAs with CIRIdeep or CIRIdeep(A) or train your own model.

## Predict

**Prediction with CIRIdeep using total RNA-seq data**

CIRIdeep provides probability of given circRNAs being differentially spliced between any of two samples. When predict with CIRIdeep, expression value of 1499 RBPs (listed in `./demo/RBPmax_totalRNA.tsv`) and splicing amount (derived from SAM alignment files) in both samples are needed. The order of RBP expression of each sample should keep exactly the same with `RBP max value file`. We recommend to process raw total RNA-seq fastq files with `CIRIquant`, which provides junction ratio of each circRNA and expression value of each gene in a one-stop manual. SAM files generated with BWA is recommended when producing splicing amount values.

```
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax_totalRNA.tsv -seqFeature ./demo/cisfeature.tsv -splicing_max ./demo/splicingamountmax_max.tsv -predict_list ./demo/predict_list.txt -model_path ./model/CIRIdeep.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_total -splicing_dir ./demo/splicingamount
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
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax_polyA.tsv -seqFeature ./demo/cisfeature.tsv -predict_list ./demo/predict_list.txt -model_path ./model/CIRIdeepA.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_polyA --CIRIdeepA
```
`--CIRIdeepA` When predict using CIRIdeepA, this parameter is needed.

Basically, the input files are similar to CIRIdeep, excluding splicing amount related files. **Notably**, the `RBP max value file` file is different from that used in CIRIdeep and all the expression values should be derived from poly(A) selected RNA-seq data. Still, when using CIRIdeep(A), the order of RBP expression of each sample should keep exactly the same with `RBP max value file`.

## Train

**CIRIdeep training**

```
python $script train -geneExp_absmax /path/to/file -seqFeature /path/to/file -splicing_max /path/to/file -outdir /out/path -RBP_dir /RBP/path -splicing_dir /splicing/path
```
Hyperparameters are given in `config.py`. `config.py` must be under the same directory with `CIRIdeep.py`. Resources are waiting to be loaded...

**CIRIdeep(A) training**

```
python $script train -geneExp_absmax /path/to/file -seqFeature /path/to/file -outdir /out/path -RBP_dir /RBP/path --CIRIdeepA
```

## Contact
Zihan Zhou. zhouzihan2018m@big.ac.cn

Please open an issue if you find bugs.

