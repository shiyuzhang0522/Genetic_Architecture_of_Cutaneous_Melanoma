#!/usr/bin/env bash
############################################################
## Rare-variant association analysis with SAIGE on UKB-RAP
##
## This script contains two modules:
##   Step 1) Fit the SAIGE null model
##   Step 2) Run chromosome-specific set-based association tests
##
## Usage:
##   export DX_PROJECT="project-xxxxxxxxxxxxxxxxxxxxxxxx"
##   bash saige_ukbrap.sh step1
##   bash saige_ukbrap.sh step2 <chr>
##
## Author: Shelley
## references: https://saigegit.github.io/SAIGE-doc/
############################################################

set -euo pipefail

: "${DX_PROJECT:?Set DX_PROJECT first (e.g., export DX_PROJECT='project-xxxx')}"

############################################################
## Step 0) Parse mode
############################################################
usage() {
  cat <<EOF
Usage:
  bash $0 step1
  bash $0 step2 <chr>

Examples:
  export DX_PROJECT="project-xxxxxxxxxxxxxxxxxxxxxxxx"
  bash $0 step1
  bash $0 step2 1
EOF
  exit 1
}

[[ $# -ge 1 ]] || usage
mode="$1"

############################################################
## Step 1) Fit SAIGE null model
############################################################
run_step1() {
  local cmd
  cmd="step1_fitNULLGLMM.R \
    --plinkFile=PLINK_for_vr_melanoma_whole_cohort \
    --phenoFile=Input_phenotype_for_nullmodel_all_melanoma_cohort.WGS.version.txt \
    --phenoCol=melanoma_status \
    --traitType=binary \
    --covarColList=PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,age_at_recruitment,genetic_sex,age_squared,sex_age,sex_age_squared,wgs_batch \
    --qCovarColList=genetic_sex,wgs_batch \
    --sampleIDColinphenoFile=eid \
    --nThreads=4 \
    --LOCO=FALSE \
    --outputPrefix=WGS_RV_melanoma_whole_cohort_null_model \
    --sparseGRMFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --isCateVarianceRatio=TRUE \
    --useSparseGRMtoFitNULL=TRUE \
    --useSparseGRMforVarRatio=TRUE \
    --SampleIDIncludeFile=WGS_shared_sample_ids.txt"

  local bedFile="${DX_PROJECT}:/2.RV_collapsing/Step1_fitnull_model/PLINK_for_vr_melanoma_whole_cohort.bed"
  local bimFile="${DX_PROJECT}:/2.RV_collapsing/Step1_fitnull_model/PLINK_for_vr_melanoma_whole_cohort.bim"
  local famFile="${DX_PROJECT}:/2.RV_collapsing/Step1_fitnull_model/PLINK_for_vr_melanoma_whole_cohort.fam"

  local subsampleFile="${DX_PROJECT}:/3.WGS_noncoding/WGS_shared_sample_ids.txt"

  local sparseGRM_matrix="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
  local sparseGRM_sample_file="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"

  local phenotype_file="${DX_PROJECT}:/3.WGS_noncoding/Input_phenotype_for_nullmodel_all_melanoma_cohort.WGS.version.txt"

  dx run app-swiss-army-knife \
    -iin="${bedFile}" \
    -iin="${bimFile}" \
    -iin="${famFile}" \
    -iin="${sparseGRM_matrix}" \
    -iin="${sparseGRM_sample_file}" \
    -iin="${phenotype_file}" \
    -iin="${subsampleFile}" \
    -icmd="${cmd}" \
    -iimage_file="${DX_PROJECT}:/docker_images/saige-latest-new-1.4.0.tar.gz" \
    --folder="${DX_PROJECT}:/3.WGS_noncoding/Step1_fitnull_model/" \
    --instance-type=mem1_ssd2_v2_x4 \
    --name=WGS_step1_fitnull_model \
    --priority=high \
    -y
}

############################################################
## Step 2) Run chromosome-specific set-based association test
## Example: melanocyte-specific cCREs
############################################################
run_step2() {
  [[ $# -eq 1 ]] || usage
  local chr="$1"

  local cmd
  cmd="step2_SPAtests.R \
    --bgenFile=/mnt/project/3.WGS_noncoding/WGS_bgen.1.2/WGS.chr${chr}.bgen.1.2version.bgen \
    --bgenFileIndex=/mnt/project/3.WGS_noncoding/WGS_bgen.1.2/WGS.chr${chr}.bgen.1.2version.bgen.bgi \
    --sampleFile=/mnt/project/3.WGS_noncoding/WGS_bgen.1.2/WGS.chr${chr}.bgen.1.2version.sample \
    --AlleleOrder=ref-first \
    --SAIGEOutputFile=noncoding.RV.melanocyte.cCRE.chr${chr}.txt \
    --chrom=${chr} \
    --LOCO=FALSE \
    --GMMATmodelFile=WGS_RV_melanoma_whole_cohort_null_model.rda \
    --varianceRatioFile=WGS_RV_melanoma_whole_cohort_null_model.varianceRatio.txt \
    --sparseGRMFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=combined_mask.mel.cCRE.chr${chr}.txt \
    --annotation_in_groupTest=CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS,GERP:GERP_JARVIS:CADD_GERP:CADD_GERP_JARVIS,JARVIS:CADD_JARVIS:GERP_JARVIS:CADD_GERP_JARVIS,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:JARVIS:GERP_JARVIS,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:GERP:GERP_JARVIS,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:JARVIS:GERP_JARVIS:GERP,CADD:CADD_GERP:CADD_JARVIS:CADD_GERP_JARVIS:JARVIS:GERP_JARVIS:GERP:Null \
    --maxMAF_in_groupTest=0.0001,0.001,0.01 \
    --is_output_markerList_in_groupTest=TRUE \
    --subSampleFile=WGS_shared_sample_ids.txt \
    --is_single_in_groupTest=TRUE \
    --is_output_moreDetails=TRUE"

  echo "${cmd}"

  local subsampleFile="${DX_PROJECT}:/3.WGS_noncoding/WGS_shared_sample_ids.txt"
  local step1_model_file="${DX_PROJECT}:/3.WGS_noncoding/Step1_fitnull_model/WGS_RV_melanoma_whole_cohort_null_model.rda"
  local step1_varianceratio="${DX_PROJECT}:/3.WGS_noncoding/Step1_fitnull_model/WGS_RV_melanoma_whole_cohort_null_model.varianceRatio.txt"
  local sparseGRM_matrix="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
  local sparseGRM_sample_file="${DX_PROJECT}:/1.GWAS-imputation/Step0.SparseGRM/sparseGRM/saige_gene_step0_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
  local target_groupfile="${DX_PROJECT}:/3.WGS_noncoding/melanocyte_cCRE_masks/combined_mask.mel.cCRE.chr${chr}.txt"
  local output_folder="${DX_PROJECT}:/3.WGS_noncoding/SAIGE.Results/melanocyte.cCRE/chr${chr}"

  dx run app-swiss-army-knife \
    -iin="${subsampleFile}" \
    -iin="${step1_model_file}" \
    -iin="${step1_varianceratio}" \
    -iin="${sparseGRM_matrix}" \
    -iin="${sparseGRM_sample_file}" \
    -iin="${target_groupfile}" \
    -icmd="${cmd}" \
    -iimage_file="${DX_PROJECT}:/docker_images/saige-latest-new-1.4.0.tar.gz" \
    --folder="${output_folder}" \
    --instance-type=mem2_ssd2_v2_x4 \
    --name="melanoma.cCRE.chr${chr}" \
    --priority=high \
    -y \
    --allow-ssh
}

############################################################
## Step 3) Dispatch
############################################################
case "${mode}" in
  step1)
    [[ $# -eq 1 ]] || usage
    run_step1
    ;;
  step2)
    [[ $# -eq 2 ]] || usage
    run_step2 "$2"
    ;;
  *)
    usage
    ;;
esac