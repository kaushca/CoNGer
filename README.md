# CoNGer
Copy Number detection in Germline samples

CoNGer is a software written mainly in Python (2.7+) and R (3.2+) to detect germline copy number variations in small-targeted capture data. Mainly, this was developed to work with data captured through Haloplex capture technology. Haloplex capture technology is different to popular hybridization capture technology, where it captures the targets using specially designed amplicons. 

CoNGer works with a large pool of germline samples. We implemented adaptive control creation, where groups of samples are selected to create pooled normal samples, in order to reduce the technological and batch based biases in data.

## General workflow

1.	Binned counts normalisation

2.	Adaptive control creation

3.	CNV segments detection and filtering

4.	Target-region CNV detection

## User guide

https://cnv.gitbook.io/conger/
