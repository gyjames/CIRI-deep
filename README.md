# CIRIdeep
CIRIdeep is a deep-learning model used to predict differentially spliced circRNAs between two biological samples using totalRNA sequencing data. An adapted version of CIRIdeep, CIRIdeepA, was trained for poly(A) selected RNA-seq data.

# Usage
The main program `CIRIdeep.py` can be used to predict differentially spliced circRNAs with CIRIdeep or CIRIdeep(A) or train your own model.

## Predict
CIRIdeep 
```
python CIRIdeep.py predict -geneExp_absmax ./demo/RBPmax.tsv -seqFeature ./demo/cisFeature.tsv -splicing_max ./demo/splicingamountmax.tsv -predict_list ./demo/sample.txt -model_path ./model/CIRIdeep.h5 -outdir ./outdir -RBP_dir ./demo/RBPexp_total -splicing_dir ./demo/splicingamount
```
