# CIRIdeep
- CIRIdeep is a deep-learning model used to predict differentially spliced circRNAs between two biological samples using totalRNA sequencing data. 
- An adapted version of CIRIdeep, CIRIdeepA, was trained for poly(A) selected RNA-seq data.

# Usage
The main program `CIRIdeep.py` can be used to predict differentially spliced circRNAs with CIRIdeep or CIRIdeep(A) or train your own model.

## Predict

**Prediction with CIRIdeep using total RNA-seq data**

CIRIdeep provides probability of given circRNAs being differentially spliced between any of two samples. When predict with CIRIdeep, expression value of 1499 RBPs (listed in `./demo/RBPmax.tsv`) and splicing amount (derived from SAM alignment files) in both samples are needed. We recommend to process raw total RNA-seq raw fastq files with `CIRIquant`, which provides junction ratio of each circRNA and expression value of each gene in a one-stop manual. SAM files generated with BWA is recommended when producing splicing amount values.

```
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax.tsv -seqFeature ./demo/cisFeature.tsv -splicing_max ./demo/splicingamountmax.tsv -predict_list ./demo/sample.txt -model_path ./model/CIRIdeep.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_total -splicing_dir ./demo/splicingamount
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

`-splicing_dir` Directory containing the splicing amount of circRNAs to be predicted in each sample. We have provided a basic script `splicing_amount.py` to produce splicing amount in samples.

**Prediction with CIRIdeep(A) using poly(A) selected RNA-seq data**

CIRIdeep(A) gives three probabilities indicating the circRNA being unchanged, having higher junction ratio in sample A or having higher junction ratio in sample B, which sum to one.
As in some cases, like in scRNA-seq or spatial transcriptomics data, only gene expression matrix is provided, splicing amount is not needed in CIRIdeep(A) any more.

```
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax.tsv -seqFeature ./demo/cisFeature.tsv -predict_list ./demo/sample.txt -model_path ./model/CIRIdeepA.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_total --CIRIdeepA
```
`--CIRIdeepA` When predict using CIRIdeepA, this parameter is needed.

Basically, the input files are similar to CIRIdeep, excluding splicing amount related files. **Notably**, the `RBPmax` file is different from that used in CIRIdeep and all the expression values should be derived from poly(A) selected RNA-seq data.

