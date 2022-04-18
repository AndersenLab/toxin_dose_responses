#!/bin/bash

bcftools view -s N2,CB4856,RC301,ECA36,ECA396,ECA248,XZ1516,MY16 /projects/b1059/analysis/WI-20200815/isotype_only/WI.20200815.hard-filter.ref_strain.vcf.gz |\\ 
    bcftools filter -i N_MISSING=0 -Oz -o Toxin_Strain_VCF.vcf.gz

tabix -p vcf Toxin_Strain_VCF.vcf.gz

plink --vcf Toxin_Strain_VCF.vcf.gz \\
            --snps-only \\
            --biallelic-only \\
            --maf 0.05 \\
            --set-missing-var-ids @:# \\
            --indep-pairwise 50 10 0.8 \\
            --geno \\
            --allow-extra-chr

        awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
        sort -k1,1d -k2,2n > markers.txt

        bcftools query -l Toxin_Strain_VCF.vcf.gz |\\
        sort > sorted_samples.txt 

        bcftools view -v snps \\
        -S sorted_samples.txt \\
        -R markers.txt \\
        Toxin_Strain_VCF.vcf.gz |\\
        bcftools query --print-header -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' |\\
        sed 's/[[# 0-9]*]//g' |\\
        sed 's/:GT//g' |\\
        sed 's/0|0/-1/g' |\\
        sed 's/1|1/1/g' |\\
        sed 's/0|1/NA/g' |\\
        sed 's/1|0/NA/g' |\\
        sed 's/.|./NA/g'  |\\
        sed 's/0\\/0/-1/g' |\\
        sed 's/1\\/1/1/g'  |\\
        sed 's/0\\/1/NA/g' |\\
        sed 's/1\\/0/NA/g' |\\
        sed 's/.\\/./NA/g' > toxin_8strain_genotype_matrix.tsv