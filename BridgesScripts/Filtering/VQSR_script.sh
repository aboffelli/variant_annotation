#!/bin/bash
# ExcessHet hard filtering
mkdir ExcessHet
for file in /home/ar7343bo-s/GATKv3.8_HC/Data/*.vcf; do gatk --java-options "-Xmx3g -Xms3g" VariantFiltration -V $file --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet -O ExcessHet/excesshet_${file##*/} 2>>ExcessHet/excesshet_error; done


# Create sites only
mkdir SitesOnly
for file in ExcessHet/excesshet*.vcf; do gatk MakeSitesOnlyVcf -I $file -O SitesOnly/sitesonly_${file#*_} 2>>SitesOnly/sites_error; done


# Training and calculate VQSLOD for SNPs
mkdir RecalTranches
#for file in SitesOnly/sitesonly*.vcf; do fname=$(echo ${file##*/}); gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator -V $file --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/ar7343bo-s/VcfTest/Resources/hapmap_3.3.b37.vcf --resource:omni,known=false,training=true,truth=true,prior=12.0 /home/ar7343bo-s/VcfTest/Resources/1000G_omni2.5.b37.vcf --resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/ar7343bo-s/VcfTest/Resources/1000G_phase1.snps.high_confidence.b37.vcf --resource:dbsnp,known=true,training=false,truth=false,prior=7 /home/ar7343bo-s/VcfTest/Resources/dbsnp_138.b37.vcf -O RecalTranches/snps_${fname%.vcf}.recal --tranches-file RecalTranches/snps_${fname%.vcf}.tranches 2>>RecalTranches/training_error; done
# using max-gaussians 6


# If generates empty files try
for file in SitesOnly/sitesonly*.vcf; do fname=$(echo ${file##*/}); 
	gatk --java-options "-Xmx3g -Xms3g" VariantRecalibrator \
	-V $file \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
	-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR \
	-mode SNP \
	--max-gaussians 1 \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /home/ar7343bo-s/Resources/hapmap_3.3.b37.vcf \
	--resource:omni,known=false,training=true,truth=true,prior=12.0 /home/ar7343bo-s/Resources/1000G_omni2.5.b37.vcf \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 /home/ar7343bo-s/Resources/1000G_phase1.snps.high_confidence.b37.vcf \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2 /home/ar7343bo-s/Resources/dbsnp_138.b37.vcf \
	-O RecalTranches/snps_${fname%.vcf}.recal \
	--tranches-file RecalTranches/snps_${fname%.vcf}.tranches 2>>RecalTranches/training_error; 
	done


# Add training and calculation for Indels
for file in SitesOnly/sitesonly*.vcf; do fname=$(echo ${file##*/}); 
	gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator \
	-V $file \
	--trust-all-polymorphic \
	-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
	-an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
	-mode INDEL \
	--max-gaussians 1 \
	-resource:mills,known=false,training=true,truth=true,prior=12.0 /home/ar7343bo-s/Resources/Mills_and_1000G_gold_standard.indels.b37.vcf \
	 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /home/ar7343bo-s/Resources/dbsnp_138.b37.vcf \
	-O RecalTranches/indels_${fname%.vcf}.recal \
	--tranches-file RecalTranches/indels_${fname%.vcf}.tranches 2>>RecalTranches/training_error; 
done

# Filter SNPs
mkdir Out
for file in ExcessHet/excesshet_*.vcf; do fname=$(echo ${file#*het_}); gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR -V $file --recal-file RecalTranches/snps_sitesonly_${fname%.vcf}.recal --tranches-file RecalTranches/snps_sitesonly_${fname%.vcf}.tranches --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode SNP -O Out/snp_vqrs_${fname} 2>>Out/filtering_error; done

# Filter Indels
mkdir Filtered
for file in Out/snp_*.vcf; do fname=$(echo ${file#*vqrs_}); gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR -V $file --recal-file RecalTranches/indels_sitesonly_${fname%.vcf}.recal --tranches-file RecalTranches/indels_sitesonly_${fname%.vcf}.tranches --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode INDEL -O Filtered/filtered_vqrs_${fname} 2>>Filtered/indel_filt_error; done
