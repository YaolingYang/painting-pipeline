# painting-pipeline

Below we described two painting pipelines.
-   Pipeline to compute haplotype components (HCs)
-   Pipeline to paint a target dataset with a reference panel

# Pipeline to compute haplotype components (HCs)
Here we describe the pipeline for computing haplotype components (HCs) through [PBWTpaint](https://github.com/richarddurbin/pbwt). We should have a genotype file from each chromosome, i.e. 22 genotype files in total. Let i denote the chromosome index, the input files are named as ``chr${i}_UKBall.vcf.gz``.

## Step 1: Run PBWTpaint for each chromosome  
The first step is to run the below command for each chromosome (i from 1 to 22, you may split it into 22 array jobs):  

``pbwt -readVcfGT chr${i}_UKBall.vcf.gz -paintSparse chr${i}_UKBall 100 2 500``

(for other input formats please just follow the pbwt instructions on reading data)  

The last parameter controls the sparsity.  

This command generates ``chr${i}_UKBall.chunklengths.s.out.gz`` for i from 1 to 22.  

## Step 2: Generate the overall chunk length matrix

The next step is to do a "weighted" sum of each entry of the N*N sparse matrix (N is the number of individuals), because PBWTpaint reports the chunk length based on the number of SNPs while it should be based on the genetic distance. Here we provide an example **C++** code ``combine_chunklength.cpp``. Please **replace row 197-198 with the number of SNPs for each chromosome in the dataset**, (row 199-201 is the genetic distance in centiMorgan for each chromosome, which do not need to change if you use a standard recombination map), and **in line 203 replace 487409 with the actual N**.  

Please ensure ``gzstream.C`` and ``gzstream.h`` are in the same directory as ``combine_chunklength.cpp``, and then compile with:

``g++ combine_chunklength.cpp -o combine -lz -lpthread -llapack -lblas -std=c++0x -g -O3``

After which you run with

``./combine``

This may take few hours to run, and afterwards we obtain the output ``full_chunklength_UKBall.txt.gz``. The first two columns are the individual indices (IND1 and IND2), and the third column is the expected chunk length that IND1 is copied from IND2 in centiMorgan. This output file represents the sparse chunk length matrix from all-vs-all painting.


## Step 3: Compute HCs via SVD

The last step is to do SVD to obtain HCs. Here we provide an example code in R, please **replace nsnp=487409 with the actual N**. 

Note that if ``full_chunklength_UKBall.txt.gz`` has more than 2^31 rows (the limit of R), then we need to remove some weakly associated individual pairs (i.e. removing the rows with the smallest numbers of the last column of ``full_chunklength_UKBall.txt.gz``) and re-weight (each row of the sparse matrix should sum up to be the total genetic distance in centiMorgan, which is 3545.04 in the standard genetic map that we use). 

``
library(data.table)
library(Matrix)
library(sparsesvd)
nsnp=487409
number_of_HCs=100
cat("Begin reading data\n")
cl <- fread("full_chunklength_UKBall.txt.gz")
cat("Begin making sparse matrix\n")
A <- sparseMatrix(cl$V1, cl$V2, x = log10(cl$V3+1), dims=c(nsnp,nsnp))
cat("Begin svd \n")
res<-sparsesvd(A,rank=number_of_HCs)
cat("Begin calculate HCs\n")
HCs <- res$u %*% diag(sqrt(res$d))
cat("Begin writing HCs\n")
colnames(HCs)=paste0(“HC”,1:number_of_HCs)
fwrite(HCs, "HCs_UKBall.csv", row.names = FALSE, sep = ',')
``

Now we get the top 100 HCs in file ``HCs_UKBall.csv``. The individual orders are the same as the original vcf.gz file.

Please also be aware that there should exist more efficient ways to do Step 2 and 3.

# Pipeline to paint a target dataset with a reference panel

Here we describe the pipeline for painting a bio-bank scale target dataset using a reference dataset on a single chromosome through [SparsePainter](https://github.com/YaolingYang/SparsePainter). Below are the files we have before painting:  

-   ``chr1_ref.vcf.gz``: The phased, non-missing reference dataset.

-   ``chr1_target.vcf.gz``: The phased, non-missing target dataset. Assume this dataset has 500,000 individuals and 50,000 SNPs.

-   ``chr1_map.txt``: The genetic map file as required by SparsePainter.

-   ``popnames.txt``: The population file of reference individuals as required by SparsePainter.

-   ``names.txt``: The names of all the target individuals.

## Step1: Split the input file into subfiles

To paint huge target samples, it is suggested to split the huge target file into small target subfiles, then we can submit multiple small jobs to run efficiently with HPC. It is usually more efficient to split phase files than vcf files, so we first convert vcf to phase files. All the commands below, unless specifically stated, are run on **Linux bash shell**.

``
pbwt -readVcfGT chr1_ref.vcf.gz -writePhase chr1_ref.phase
pbwt -readVcfGT chr1_target.vcf.gz -writePhase chr1_target.phase
``

Then we create a folder ‘chr1split’ to store the subfiles.

``mkdir chr1split``

Now we run below commands in **Python** to split target files into 500 subfiles, each subfile contains 1000 individuals. Note that this is an example code, more efficient codes for splitting files are possible.

```
import os
chr = '1'
input_file = f"chr{chr}_target.phase"
total_files = 500
lines_per_file = 2000

def split_file(input_file, total_files, lines_per_file, chr):
    with open(input_file, 'r') as file:
       next(file)  # Skip the first line
       second_line = next(file)
       third_line = next(file)

    with open(input_file, 'r') as file:
        for _ in range(3):
            next(file)  

        file_number = 1
        line_count = 0
        output_file = gzip.open(f"chr{chr}split/chr{chr}_target{file_number}.phase", 'wt')
        output_file.write(f"{lines_per_file}\n")  
        output_file.write(second_line)           
        output_file.write(third_line)             

        for line in file:
            if line_count == lines_per_file:
                output_file.close()
                file_number += 1
                line_count = 0
                output_file = gzip.open(f"chr{chr}split/chr{chr}_target{file_number}.phase", 'wt')
                output_file.write(f"{lines_per_file}\n")  
                output_file.write(second_line)
                output_file.write(third_line)

            output_file.write(line)
            line_count += 1

        output_file.close()

split_file(input_file, total_files, lines_per_file, chr)
```

To save storage space, we could compress the reference file, and remove ``chr1_target.phase`` which has already been split into subfiles.

``
gzip chr1_ref.phase
rm chr1_target.phase
``

Also, we should update the name files for the individuals in each subfile. These files are stored in the folder ``namefile``.  

```
input_file="names.txt"
output_dir="namefile"
mkdir -p "$output_dir"

lines_per_file=1000
file_count=1
line_count=0

while IFS= read -r line; do
    if [ $line_count -eq 0 ]; then
        output_file="${output_dir}/target${file_count}.txt"
        > "$output_file"
    fi

    echo "$line" >> "$output_file"

    ((line_count++))

    if [ $line_count -eq $lines_per_file ]; then
        ((file_count++))
        line_count=0
        echo $file_count
    fi
done < "$input_file"
```

Now we have finished data preprocessing. ``chr1_target.vcf.gz`` has been split into 1000 subfiles in phase format in folder ``chr1split``; ``names.txt`` has been split into 1000 subfiles in folder ``namefile``. ``chr1_ref.vcf.gz`` has been converted to ``chr1_ref.phase.gz``. Then the painting starts.

## Step2: Generate weight file through ref-vs-ref painting

Here we generate the weight file for each reference sample at each SNP. We first generate refnames.txt which is the names for reference individuals:

``
cut -d ' ' -f1 popnames.txt > refnames.txt
``

Then we do ref-vs-ref painting (without leave-one-out) to generate the weight file:

``
 ./SparsePainter -weight -reffile chr1_ref.phase.gz -targetfile chr1_ref.phase.gz -popfile popnames.txt -mapfile chr1_map.txt -namefile refnames.txt -indfrac 1 -out chr1_ref
``

This command generates ``chr1_ref_weight.txt.gz``, which is the weight file, and ``chr1_ref_fixedlambda.txt``, which is the estimated recombination scaling constant. If you paint multiple chromosomes, and assume the estimated recombination scaling constant is 100, please add command ``-fixlambda 100`` in the above command when generating the weight file for the other chromosomes.


Step3: Perform the corrected painting for all the target individuals

Before we paint all the target individuals, we need to determine the recombination scaling constant lambda, and use the fixed lambda to paint all the individuals. To estimate lambda, we only need to paint one target subset:

``
mkdir chr1
 ./SparsePainter -reffile chr1_ref.phase.gz -targetfile chr1split/chr1_target1.phase -weightfile chr1_ref_weight.txt.gz -popfile popnames.txt -mapfile chr1_map.txt -namefile namefile/target1.txt -indfrac 1 -prob -chunklength -chunkcount -probstore linear -out chr1/chr1_target1
``

This generates ``chr1/chr1_target1_fixlambda.txt`` and other files (described below). Assume the estimated recombination scaling constant is 100 from ``chr1/chr1_target1_fixlambda.txt``, then we use this fixed lambda to paint all the other subfiles (and other chromosomes if it applies). We usually submit the remaining 499 array jobs on HPC. Let ``SLURM_ARRAY_TASK_ID`` denote the array task IDs, then we run the below command to paint all the subfiles:

``
 ./SparsePainter -reffile chr1_ref.phase.gz -targetfile chr1split/chr1_target${SLURM_ARRAY_TASK_ID}.phase -popfile popnames.txt -weightfile chr1_ref_weight.txt.gz -mapfile chr1_map.txt -namefile namefile/target${SLURM_ARRAY_TASK_ID}.txt -fixlambda 100 -prob -chunklength -chunkcount -probstore linear -out chr1/chr1_target${SLURM_ARRAY_TASK_ID}
``

It is optional to paint with ``-LDAS -AAS -aveSNP -aveind``, etc.

Now we finish the painting, and get the below output files for each target subfile indexed ``SLURM_ARRAY_TASK_ID``.

``chr1/chr1_target${SLURM_ARRAY_TASK_ID}_prob.txt.gz``: The local ancestry probabilities for each target sample at each SNP stored in linear form (see SparsePainter manual).
``chr1/chr1_target${SLURM_ARRAY_TASK_ID}_chunklength.txt.gz``: The expected length (in centiMorgan) of copied chunks of each local ancestry for each target sample.
``chr1/chr1_target${SLURM_ARRAY_TASK_ID}_chunkcount.txt.gz``: The expected number of copied chunks of each local ancestry for each target sample.

## Step4: Post-processing output files

Please follow the [instructions](https://github.com/YaolingYang/SparsePainter) to merge these subfiles, extract paintings for analysis, etc. Below we provide an example in **Python** to merge the probability (``chr1_target_prob.txt.gz``), chunk length (``chr1_target_chunklength.txt.gz``) and chunk count (``chr1_target_chunkcount.txt.gz``) files.

```
import gzip
import os
import sys

chr = '1'
total_files = 500

#Function to merge files
def merge_probfiles(file_type, output_file_name):
    with gzip.open(output_file_name, 'wt') as output_file:
        for i in range(1, total_files+1):
            print(f"Processing file {i} for {file_type}")
            full_path = f"chr{chr}/chr{chr}_target{i}_{file_type}.txt.gz"

            if os.path.isfile(full_path):
                with gzip.open(full_path, 'rt') as input_file:
                    if i != 1:
                        next(input_file, None)
                        next(input_file, None)
                    for line in input_file:
                        output_file.write(line)
            else:
                print(f"File not found: {full_path}")

def merge_chunkfiles(file_type, output_file_name):
    with gzip.open(output_file_name, 'wt') as output_file:
        for i in range(1, total_files+1):
            print(f"Processing file {i} for {file_type}")
            full_path = f"chr{chr}/chr{chr}_target{i}_{file_type}.txt.gz"

            if os.path.isfile(full_path):
                with gzip.open(full_path, 'rt') as input_file:
                    if i != 1:
                        next(input_file, None)

                    for line in input_file:
                        output_file.write(line)
            else:
                print(f"File not found: {full_path}")

#Merge 'prob' files
merge_probfiles('prob', f"chr{chr}_target_prob.txt.gz")

#Merge 'chunklength' files
merge_chunkfiles('chunklength', f"chr{chr}_target_chunklength.txt.gz")

#Merge 'chunkcount' files
merge_chunkfiles('chunkcount', f"chr{chr}_target_chunkcount.txt.gz")
```
