## Example script notes:
#  1)  Samples in the example were downsampled and only contain chromosome 17, in order to make the file small enough to follow github restriction
#  2)  This example is primarily intended for verifying installation, and demonstrated required inputs and content formats  
#  3)  Reference assembly (FASTA file) for the provided example can be obtained from IGSR website, which based on GRCh37  
#			http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
#
#  Running example:
#     - assumes the working directory is at testrun, and outputting to the current folder
#     - assumes reference fasta file is in the working directory, otherwise adjust the path accordingly
#
#  1) Control pool creation (outputs to folder named ctrl):
#		python ../../create_pool.py -c controls.txt -b ../testdata/example_Halo_2_annotated_chr17.bed --out ctrl --fasta human_g1k_v37.fasta --info controls-info.tsv
#
#  2) CNV analysis (outputs to folder named tumour-cnv):
#		python ../../CoNGer.py -t tumour.txt -c ctrl -b ../testdata/example_Halo_2_annotated_chr17.bed --out tumour_cnv --fasta human_g1k_v37.fasta --info tumour-info.tsv
       

