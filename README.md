# DAPseq_TFanalysis
This is the compagnion repository site for book chapter YYYY.

The main scripts pipeline.sh analyze DAP-seq sequence reads to call peaks and build and evaluate transcription factors binding site model.

All pipeline dependencies but MEME and MSPC are packed into the DAPseqEASY conda environment.

Setting up the DAPseqEASY conda environment

conda env create -f DAPseqEASY.yml

pipeline dependencies not included in the DAPseqEASY conda environment that you'll have to install locally.

meme version 4.12.0

MSPC version 5.4.0

from you working directory

mkdir scripts

move  pipeline.sh, meme2pfm.sh, ROC.R and scores.py into the scripts directory

conda activate DAPseqEASY

./pipeline.sh
