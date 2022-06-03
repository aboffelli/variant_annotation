# README file

GitHub repository: https://github.com/aboffelli/variant_annotation

## Setting up the softwares and environments

### Installing GATK

After downloading the file, and unzipping it, there is no installation needed to use the program. To facilitate the command line we can add the directory to our PATH.

```shell=
## Download and unzip the program
wget https://github.com/broadinstitute/gatk/releases/download/4.2.2.0/gatk-4.2.2.0.zip
unzip gatk-4.2.2.0.zip

## All necessary files to run the program will be inside the gatk-4.2.2.0 directory. So we can add this directory to the path. By adding the following line to the .bashrc file. Remember to change the path according to where your directory is saved.
export PATH=$PATH:~/bin/gatk-4.2.2.0/

```

### VCFTools and VEP
We can install both software with conda. The following line creates a new conda environment called ***e_vep***.
```shell=
## Create the new conda evironment with all necessary packages.
conda create --name e_vep -c bioconda ensembl-vep mysql perl-bio-bigfile perl-bio-db-hts vcftools perl-vcftools-vcf

## Activate the new environment
conda activate e_vep
```

We also need to download the ensembl-vep package to be able to download the cache.

```shell=
## Download the package
git clone https://github.com/Ensembl/ensembl-vep.git

## Use the instalation script
cd ensembl-vep
perl INSTALL.pl
## Follow the instructions on the screen and download the cache for Homo sapiens.
```

### Installing MMSplice

We install MMSplice in a different environment.
```shell=
## Create and activate the new environment
conda create -n splice
conda activate splice

## download the MMSplice package
git clone https://github.com/gagneurlab/MMSplice

cd MMSplice/VEP_plugin
cp MMSplice ~/.vep/Plugins/

## Install Cython and mmsplice
pip install Cython
pip install mmsplice

## update tensorflow, numpy and scipy
pip install --upgrade tensorflow
pip install --upgrade numpy
pip install --upgrade scipy

```

### Installing PLINK
Since plink uses Python 2 we install in a different environment

```shell=
## Create a new environment and install plink
conda create -n plink_env plink
```

### Resource files for VEP

All resource files can be stored in a directory.
```shell=
mkdir Resources

cd Resources
```
#### PhyloP and GERP resource files
```shell=
## PhyloP
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/phyloP100way/hg19.100way.phyloP100way.bw

## GERP
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/All_hg19_RS.bw.
```

The Swegen resource is not publicly available.

#### Encode resouce for RBP binding sites.
Download the file list
```shell=
mkdir EncodeFiles
cd EncodeFiles

wget https://www.encodeproject.org/batch_download/?control_type!=*&status=released&perturbed=false&assay_slims=RNA+binding&assembly=hg19&files.file_type=bed+narrowPeak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=eCLIP&type=Experiment
```

The first line in files.txt is the metadata link, which can be downloaded
```shell=
wget https://www.encodeproject.org/metadata/?control_type%21=%2A&status=released&perturbed=false&assay_slims=RNA+binding&assembly=hg19&files.file_type=bed+narrowPeak&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=eCLIP&type=Experiment
```

After downloading the file list and the metadata, select only the files with assembly hg19 and Biological replicate 1,2.

```sh
grep "GRCh37" ENCODE_RBP_metadata.tsv | grep "1, 2" > files_to_download.tsv 
```

The first column in the metadata is the file accession, use that to download only the needed files.

```sh
cat ~/files_to_download.tsv | cut -f1 | while read line; do file=$(grep "$line" ~/ENCODE_RBP_files.txt); wget $file ; done

## Decompress the bed files.
for file in *; do filename=${file%.gz}; gunzip -c $file > $filename; done
```

Prepare the bed file for VEP

```sh
## Join only column 1-8 from all files, and remove chr from the chromosome column and _IDR from prot column.
for f in *.bed; do
     while read line; do
         feature=$(echo "$line" | cut -f4-8 | tr $'\t' ':');
         pos=$(echo "$line" | cut -f1-3);
         echo "$pos"$'\t'"$feature" | sed 's/_IDR//';
     done < "$f" >> encode_rna_binding.bed;
 done

## Sort the file and compress
cat encode_rna_binding.bed | sort -k1,1 -k2,2n -k3,3n | bgzip -c > encode_rna_binding_1.bed.gz

## Index with tabix
tabix -p bed encode_rna_binding_1.bed.gz
```

#### Clinvar resources
```shell=
## Download the tab delimited ClinVar file.
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

## Retrieve only GRCh37 lines and save in a file.
zcat variant_summary.txt.gz | grep "GRCh37" > variant_summary_GRCh37.txt
```


## SWEA analysis
We first create a main directory where we will run all scripts.
```shell=
mkdir SWEA
cd SWEA
```
### Filtering
With all the VCF files in a directory called **Data**, we can run the filtering script.

```shell=
## In the main directory where the directory Data is found use the following line.
bash hard_filtering_script.sh
```

Now we can create the QC plots for the filtering.
```shell=
bash table_script.sh
```

The resulting tables can be ploted with *filtered_var.R*.

### Annotation
Now we run all the following order of scripts in the main directory. 

#### VEP annotation
```shell=
## Activate conda env for VEP.
conda activate e_vep

## NOTE: These steps may take hours or days, run the scripts with nohup or using software that you can detach, such as tmux or screen. 

# First annotation with VEP
bash annotation_script.sh
```

We can create tables for evaluating how many variants are on target and off target
```shell=
## First we need a reference gff file for gene positions
cd ~/Resources
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh37_mapping/gencode.v38lift37.annotation.gff3.gz

## And run the python script in the main directory supplying a file with the genes sequenced in the study.
cd ~/SWEA
python variant_position_parser.py ~/Resources/gencode.v38lift37.annotation.gff3.gz whole_gene_list.txt

## And generate the tables
cd Annotation/
bash annotation_table.sh

## This script retrieves the number of variants that failed a specific filter, we used for SOR3 and QD2S.
bash filter_number_fail_script.sh 
```

After plotting with *annotation_plots.R* we realized we could change the SOR value to > 4.
```shell=
## We run the script to make the changes.
cd ~/SWEA
python sor3_to_sor4.py
```

We can repeat the plotting steps for the edited files, by running the scripts in the `Edited/` directory.

We then follow the annotation steps.

```shell=
## Encode annotation with VEP.
bash rna_binding_annotation.sh

## Gene name annotation with VEP.
bash gene_name_annotation.sh
```

#### Custom annotation
For the custom annotation, we need files containing ESE and ESS hexamers in our `Resources/`directory, the papers used are cited in the report.
```shell=
## First custom annotation ** Need the ESE ESS files.
python custom_annotation.py
```

Here we can generate more tables to plot with `count_ESS_ESE_RBP.py` for ESE/ESS and RBP, and `variant_distribution_table.py` for deltaRSCU, conservation and known and novel distribution. R plotting with scripts `es_rbp_counter.R`, `variant_distribution.R`, and `conserv_rscu_plots.R`.

Since MMSplice takes a long time to run we can put it to run in the backgroud while continue the other annotation.
```shell=
## Activate conda environment for MMSplice. 
conda activate splice
bash mmsplice_script.sh
```

#### ClinVar annotation
```shell=
## ClinVar custom python script.
python clinvar.py
```

In this step we can generate tables and plots from ClinVar annotation with `clinvar_tables.py` and`clinvar_plots.R`.

#### Adding MMSplice to ClinVar annotated files
After MMSplice finish running we can add the results to the fully annotated files
```shell=
python mmsplice_parser.py
```

#### Removing samples with pathogenic variants
Create a file with all the samples that contain pathogenic variants (the first column of the samples_pathogenic.txt table). Then in a copy of the directory containing all samples, remove the ones that are in the text file.

```sh
cat ClinVarTables/samples_pathogenic.txt | cut -f1 | sort -u > ../samples_with_patogenic.txt

## Remove the samples from the new directory with all samples. 
cat samples_with_patogenic.txt | while read line; do rm MMSamplesWithoutPathogenic/mmsplice_clinvar_custom_gene_name_encode_edited_fs_filtered_${line}_VEP.vcf; done
```

#### Synonyous table
Finally we can generate the resulting synonymous table
```shell=
python synonymous_table.py

## Also create a table for all variants
python complete_table.py

## The calculation of the allele frequency in synonymous_table can be updated for all samples using the complete table
## Change the AF SWEA in synonymous_table

cut -f1,2 synonymous_table.txt | while read line; 
do 
    old=$(grep "$line" synonymous_table.txt | cut -f7); 
    new=$(grep "$line" complete_table.txt | cut -f6); 
    grep "$line" synonymous_table.txt | sed "s/$old/$new/" >> new_synonymous_table.txt; 
done
```

#### Interactive heatmap
The interactive heatmap is created using the `new_synonymous_table.txt` in the R script `synonymous_plot.R`.


## BRIDGES
Bridges files have an original filtering that we are changing later, so we can annotate the files first.

### Annotation
```shell=
conda activate e_vep

## NOTE: These steps may take weeks, run the scripts with nohup or using software that you can detach, such as tmux or screen.

## Run the first annotation script with VEP.
bash annotation_script.sh

## Followed by the Encode annotation with VEP.
bash rna_binding_annotation.sh

## Now we create the directory structure for the next annotation.
rsync -av -f"+ */" -f"- *" EncodeAnnotation/ CustomAnnotation/

## And run the next script
python custom_annotation.py

## Now we create the directory structure for the next annotation.
rsync -av -f"+ */" -f"- *" CustomAnnotation/ ClinVar/

## And run the next script
python clinvar.py
```

### Plotting

Since we have all the annotation, we can make the tables and plot in R.
```shell=
## Run the script to generate the tables for known and novel, and ClinVar plots.
python novel_known_table.py ClinVar/

## For the ClinVar tables script we need to go to the right directory
cd ClinVar
python clinvar_tables.py
```

For plotting use the resulting files in the scripts *bridges_known_novel_plots.R* and *clinvar_plots.R*


### Additional filtering
Now that we realized that the BRIDGES files need more filtering, we need to change the already annotated files.
```shell=
## So we prepare the directory structure
rsync -av -f"+ */" -f"- *" ClinVar/ FilteredClinVar/

## And run the script.
python original_filtering.py
```

Now we can plot the filtered values by running the same commands as before, but changing the directory to FilteredClinVar.
```shell=
## Run the script to generate the tables for known and novel, and ClinVar plots.
python novel_known_table.py FilteredClinVar/

# For the ClinVar tables script we need to go to the right directory
cd FilteredClinVar
python clinvar_tables.py
```

### Association analysis
PLINK requires an specific input, so first of all we run the script to create the input files for PLINK.

The usage of plink_input.py script is described below:
```
Usage: plink_input.py [N] [--family] [--synonymous] [-h]

    N: positional argument integer with the cut-off number for the number of samples containing a variant
    --family: Boolean to remove samples that do not have family history
    --synonymous: Boolean to remove non-synonymous variants
    --pathogenic: Boolean to remove samples with pathogenic variants
    -h --help: Prints this help message and exits.
```

The file *BRIDGES_fluidigm_juno_panel_primers.txt* must be in the Resources directory for this script.
```shell=
## We create two options, one with synonymous variants of all samples without pathogenic variants. 
python plink_input.py 1 --synonymous --pathogenic

## And other with only samples with family history without pathogenic variants.
python plink_input.py 1 --synonymous --pathogenic --family
```

Now we can run PLINK.
```shell=
## Since we want to run PLINK for both files we can use a loop.
for f in *.map; do plink --file ${f%.map} --no-fid --no-parents --assoc --adjust --out plink_${f%.map}; done
```
    --no-fid: ped file without family id column
    --no-parents: ped file without parents id column
    --assoc: performs case-control association
    --adjust: outputs a file with adjusted p-values
    
We retrieve the significant synonymous variants in a table.
```shell=
python association_synonymous_table.py plink_filt1_synonymous_no_pathogenic.assoc.adjust

python association_synonymous_table.py plink_filt1_family_synonymous_no_pathogenic.assoc.adjust
```

#### Counting the overlapping variants 
```shell=
cat filt1_synonymous_no_pathogenic_table.txt filt1_family_synonymous_no_pathogenic_table.txt | grep -v "#" | cut -f1,2| uniq | sort | uniq -c | grep "2 " | wc -l
```

#### Get variants in common with SWEA.
```shell=
## Samples with family history
cat raidset/FilteredClinVar/filt1_family_synonymous_no_pathogenic_table.txt | cut -f1,2 | while read var; 
do 
    pos=$(echo ${var} | cut -d ' ' -f1); 
    base=$(echo ${var} | cut -d ' ' -f2); 
    grep -P "${pos}\t${base}" bridges_synonymous_table.txt >> variants_in_common_S_B_family.txt; 
done

## All samples.
cat raidset/FilteredClinVar/filt1_synonymous_no_pathogenic_table.txt | cut -f1,2 | while read var; 
do 
    pos=$(echo ${var} | cut -d ' ' -f1); 
    base=$(echo ${var} | cut -d ' ' -f2); 
    grep -P "${pos}\t${base}" bridges_synonymous_table.txt >> variants_in_common_S_B_all_samples.txt; 
done
```

#### Add OR in the tables
```shell=
bash or_parser.sh variants_in_common_S_B_family.txt plink_bridges_filt1_family_synonymous_no_pathogenic.assoc

bash or_parser.sh variants_in_common_S_B_all_samples.txt plink_bridges_filt1_synonymous_no_pathogenic.assoc

bash or_parser.sh filt1_synonymous_no_pathogenic_table.txt plink_bridges_filt1_synonymous_no_pathogenic.assoc

bash or_parser.sh filt1_family_synonymous_no_pathogenic_table.txt plink_bridges_filt1_family_synonymous_no_pathogenic.assoc
```

Lastly we can add the calculation of allele frequency from BRIDGES files in the synonymous table from SWEA.
```shell=
python synonymous_comparison.py
```
