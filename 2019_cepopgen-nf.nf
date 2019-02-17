#!/usr/bin/env nextflow 

date = new Date().format( 'yyyyMMdd' )


params.help                   = null
params.config                 = null
params.out_base               = null
params.anc                    = null
params.ref                    = params.reference
params.cpu                    = "4"
params.snv_vcf                = params.snv_vcf
params.sv_vcf                 = params.sv_vcf
params.admix_k                = params.admix_k
params.admix_ld               = params.admix_ld
params.pops                   = params.pops
params.out                    = "${date}-${params.out_base}"


K = Channel.from(2..params.admix_k)


ld_range = Channel.from(params.admix_ld)
// use this to test more than one ld. 
// format of admix_ld = "ld\n0.1\n0.2\n0.3"
//ld_range = Channel.from(params.admix_ld)
//  .splitCsv(header:true)
//  .map{ row -> row.ld }
// lds = params.admix_ld.replaceAll(/\n/, ",").replaceAll(/ld,/, "")

log.info ""
log.info "----------------------------------------------------------------"
log.info "        C. elegans PopGen Paper Nextflow Pipeline "
log.info "----------------------------------------------------------------"
log.info ""


if (params.help) {
    log.info "----------------------------------------------------------------"
    log.info "                      USAGE                                     "
    log.info "----------------------------------------------------------------"
    log.info ""
    log.info "nextflow cepopgen-nf.nf --out_base Analysis --anc XZ1516"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--out_base             String                Name of folder to output results"
    log.info "--anc                  String                Name of ancestor strain to use"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "Optional arguments:"
    log.info "Information describing the stucture of the input files can be located in input_files/README.txt"
    log.info ""
    log.info "--ref                  FILE                 Location of the reference genome to use"
    log.info "--config               FILE                 Location of a custom configuration file"    
    log.info "--snv_vcf              FILE                 Location to the small variant VCF to use for analysis"
    log.info "--sv_vcf               FILE                 Location to the structural variant VCF to use for analysis"
    log.info "--pops                 FILE                 Location to sample population file"    
    log.info "--admix_k              INTEGER              Maximum admixture population size to test"
    log.info "--admix_ld             STRING               Comma separated string of LDs to test for admixture analysis"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    log.info "--email                STRING               email address for job notifications"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    log.info "--------------------------------------------------------"
    log.info ""
    log.info " Required software packages to be in users path"
    log.info "BCFtools               v1.9"
    log.info "plink                  v1.9"
    log.info "RAxML-ng               v0.5.1b"
    log.info "VCFtools               v0.1.16"
    log.info "ADMIXTURE              v1.3"
    log.info "vcfanno                v0.2.8"
    log.info "EIGENSOFT              v6.1.4"
    log.info "VARISCAN               v2.0"
    log.info "--------------------------------------------------------"    
    exit 1
} else {



log.info ""
log.info "Reference                               = ${params.ref}"
log.info "Ancestor                                = ${params.anc}"
log.info "Small Variant VCF                       = ${params.snv_vcf}"
log.info "Structural Variant VCF                  = ${params.sv_vcf}"
log.info "Population File                         = ${params.pops}"
log.info "Max Population                          = ${params.admix_k}"
log.info "LD to test                              = ${params.admix_ld}"
log.info "cpu                                     = ${params.cpu}"
log.info "Output folder                           = ${params.out}"
log.info "GFF3 File                               = ${params.gff}"
log.info ""
}

/*
~ ~ ~ > * Define Contigs 
*/
CONTIG_LIST = ["I", "II", "III", "IV", "V", "X"]
Channel.from(CONTIG_LIST)
       .into{contigs_popgenome_window;
             contigs_popgenome_gene;
             contigs_h21;
             contigs_raisd;
             contigs_popgenome_gene_subpops;
             contigs_other1;
             contigs_other2}

// initialize population 
File pop_file = new File(params.pops);

pop_strains = Channel.from(pop_file.collect { it.tokenize( ' ' ) })
             .map { POP, MAF, SM -> [POP, MAF, SM] }

pop_strains
  .combine(ld_range)
  .into{subset_vcf_input;
        plink_input; 
        print_plink_input }

// initialize VCF channel and spread to different channels for analysis
small_vcf = Channel.fromPath(params.snv_vcf)
small_index = Channel.fromPath("${params.snv_vcf}" + ".tbi")


// move this below VCF annotation when we decide to switch analysis to using masked VCF.

small_vcf
  .spread(small_index)
  .into { smallvcf_ancestor;;
          smallvcf_annotations;
          }

// initialize eigenstrat parameter file chanel
eigenstrat_noremoval = Channel.fromPath(params.eigen_par_no_removal)
eigenstrat_removal = Channel.fromPath(params.eigen_par_outlier_removal)

/*
==================================
~ > *                        * < ~
~ ~ > *                    * < ~ ~
~ ~ ~ > *  ANNOTATE VCF  * < ~ ~ ~
~ ~ > *                    * < ~ ~
~ > *                        * < ~
==================================
*/


/*
------------ Extract ancestor strain from the VCF and make bed file for annotations 
*/

process extract_ancestor_bed {

    publishDir "${params.out}/ANNOTATE_VCF", mode: 'copy'

    cpus 1

    input:
      set file(vcf), file(vcfindex) from smallvcf_ancestor

    output:
      set file("ANC.bed.gz"), file("ANC.bed.gz.tbi") into anncestor_bed

      """

        bcftools query --samples ${params.anc} -f '%CHROM\\t%POS\\t%END\\t[%TGT]\\n' ${vcf} |\\
        awk -F"/" '\$1=\$1' OFS="\\t" |\\
        awk '{print \$1, \$2 = \$2 - 1, \$3, \$4}' OFS="\\t" |\\
        bgzip > ANC.bed.gz

        tabix ANC.bed.gz
        echo "ANCESTOR DONE"
      """
}

/*
------------ Annotate small variant VCF  
*/

process annotate_small_vcf {

    publishDir "${params.out}/ANNOTATE_VCF", mode: 'copy'

    cpus 1

    input:
      set file(vcf), file(vcfindex) from smallvcf_annotations
      set file("ANC.bed.gz"), file("ANC.bed.gz.tbi") from anncestor_bed

    output:
      set file("Ce330_annotated.vcf.gz"), file("Ce330_annotated.vcf.gz.tbi") into annotated_vcf

      """

        vcfanno ${params.vcfanno_config} ${vcf} |\\
        awk '\$0 ~ "#" || \$0 !~ "Masked" {print}' |\\
        bcftools filter -i N_MISSING=0 -Oz -o Ce330_annotated.vcf.gz

        tabix -p vcf Ce330_annotated.vcf.gz

        echo "done it again yes1"
      """
}

annotated_vcf
  .into { annotated_vcf_haplotype;
          annotated_vcf_popgen;
          annotated_vcf_eigenstrat;
          smallvcf_admixture;
          smallvcf_phylo;
          smallvcf_haplotype;
          smallvcf_recombination;
          smallvcf_classic_popgen_window;
          smallvcf_classic_popgen_gene;
          smallvcf_classic_popgen_gene_admix_pops;
          smallvcf_haplotype_popgen;
          smallvcf_popstructure;
          smallvcf_eigenstrat;
          smallvcf_dapc;
          smallvcf_h21;
          smallvcf_raisd;
          smallvcf_other}

/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > *  ADMIXTURE ANALYSIS  * < ~ ~ ~
~ ~ > *                          * < ~ ~ 
~ > *                              * < ~
========================================
*/

smallvcf_admixture
  .spread(subset_vcf_input)
  .into { vcf_admix_plink;
          vcf_admix_print }


/*
------------ Prune VCF, Generate PLINK files
*/

process vcf_to_ped {

    tag {"PRUNING VCF FOR ADMIXTURE"}

    publishDir "${params.out}/ADMIXTURE/PLINK/", mode: 'copy'

    input:
      set file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld) from vcf_admix_plink

    output:
      set file("ce_norm.vcf.gz"), file("ce_norm.vcf.gz.tbi"), val(nSM), val(maf), val(samples), val(ld), file("*.map"), file("*.ped") into plink_output
      file("plink.prune.in") into pruned_marker_set

      """

      bcftools norm -m +snps ${vcf} -Oz -o ce_norm.vcf.gz
      tabix -p vcf ce_norm.vcf.gz

      plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --maf ${maf} --set-missing-var-ids @:# --indep-pairwise 50 10 ${ld} --allow-extra-chr 
      
      plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --maf ${maf} --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out LD_${ld}_MAF_${maf} --allow-extra-chr 
      """

}

/*
------------ Generate random seed for running ADMIXTURE 10 times
*/

Channel.from( 1000..2000 )
       .randomSample( 10, 234 )
       .into{ rseed_find_best_k;
              rseed_run_best_k}


/*
------------ Append population range and random seed to pruned plink files for ADMIXTURE analysis
*/

plink_output
  .into{ admixture_find_best_k;
         admixture_run_best_k}

K
  .spread(rseed_find_best_k)
  .spread(admixture_find_best_k)
  .into{ admixture_genome;
         admixture_regions}

/*
------------ Run ADMIXTURE analysis
*/

process run_admixture {

    tag { " ${pop} - ${rseed} "}

    publishDir "${params.out}/ADMIXTURE/${pop}/", mode: 'copy', pattern: '*.P'
    publishDir "${params.out}/ADMIXTURE/${pop}/", mode: 'copy', pattern: '*.Q'
    publishDir "${params.out}/ADMIXTURE/${pop}/", mode: 'copy', pattern: 'log_*'

    cpus 4

    input:
      set pop, rseed, file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped) from admixture_genome

    output:
      set pop, rseed, file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped), file("log_${pop}_${rseed}.out"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.P"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.Q") into admixture_output
      set pop, rseed, file("log_${pop}_${rseed}.out"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.P"), file("LD_${ld}_MAF_${maf}_${pop}_${rseed}.Q") into admix_results

      """

        admixture --cv=10 -s ${rseed} ${ped} ${pop} -j4 | tee log_${pop}_${rseed}.out

        mv LD_${ld}_MAF_${maf}.${pop}.P LD_${ld}_MAF_${maf}_${pop}_${rseed}.P
        mv LD_${ld}_MAF_${maf}.${pop}.Q LD_${ld}_MAF_${maf}_${pop}_${rseed}.Q
      """
}

/*
------------ Combine 10 independent runs of each K
*/

admix_results
  .groupTuple()
  .into{ grouped_admix;
         grouped_print}

//grouped_print.println()

/*
------------ Combine log files of 10 independent runs of each K 
*/

process concat_replicate_logs {

    tag { " ${pop} "}

    echo true

    executor 'local'

    input:
      set pop, rseed, file(admixlog), file(admixp), file(admixq) from grouped_admix

    output:
      file("K${pop}_summary.txt") into concatenated_logs


    """
      grep -h CV log*.out |\\
      cut -f3- -d" " |\\
      sed 's/(\\|)\\|:\\|K=//g' > K${pop}_summary.txt
    """

}

/*
------------ Process CV results from ADMIXTURE analysis
*/

process concat_pop_logs {

    publishDir "${params.out}/ADMIXTURE/CV_Summary/", mode: 'copy'

    echo true

    executor 'local'

    input:
      file(clog) from concatenated_logs.toSortedList()

    output:
      set file("admix_replicates_CV.tsv"), file("admix_summarized_CV.txt"), file("bestK.txt") into cv_summary
      file("bestK.txt") into bestk

    """

      # FULL RESULTS
      cat *summary.txt |\\
      sort -k1n |\\
      awk '\$1=\$1' OFS="\\t" |\\
      awk 'BEGIN{OFS="\\t"; print "K", "CV"}; {print \$0} OFS="\\t"' > admix_replicates_CV.tsv

      # Means of replicates
      cat *summary.txt |\\
      sort -k1n |\\
      awk '\$1=\$1' OFS="\\t" |\\
      awk 'BEGIN{OFS="\\t"; print "K", "CV"}; {print \$0} OFS="\\t"' |\\
      datamash -g 1 mean 2 -H |\\
      sed 's/GroupBy(\\|)\\|mean(//g' > admix_summarized_CV.txt

      # FIND BEST K - FIRST K WHERE NEXT HIGHER K CV VALUE IS SAME TO 2 decimal places
      cat *summary.txt |\\
      sort -k1n |\\
      awk '\$1=\$1' OFS="\\t" |\\
      awk 'BEGIN{OFS="\\t"; print "K", "CV"}; {print \$0} OFS="\\t"' |\\
      datamash -g 1 mean 2 -H |\\
      sed 's/GroupBy(\\|)\\|mean(//g' |\\
      awk 'NR>1{print \$0, sprintf("%3.2f", \$2-p)} {p = \$2}' |\\
      sed 's/-//g' |\\
      awk '\$3 == 0.00 {print}' |\\
      head -1 |\\
      cut -f-1 > bestK.txt

      # ADD PLUS MINUS 2 RANGE TO BEST K
      bk=`head -1 bestK.txt`

      START=\$(( 1+bk ))
      END=\$(( 2+bk ))
      for ((i=START;i<=END;i++)); do
          echo \$i >> bestK.txt
      done

      START=\$(( bk-2 ))
      END=\$(( bk-1 ))
      for ((i=START;i<=END;i++)); do
          echo \$i >> bestK.txt
      done

    """
}

/*
------------ Re-Run ADMIXTURE Analysis using Best K value 
*/

bestk
  .splitText()
  .set{k_range}

admixture_run_best_k
  .spread(k_range)
  .into{ rerun_admixture;
         rerun_admixture_extra}

process run_admixture_besk_k {

    publishDir "${params.out}/ADMIXTURE/BEST_K/", mode: 'copy', pattern: '*.P'
    publishDir "${params.out}/ADMIXTURE/BEST_K/", mode: 'copy', pattern: '*.Q'

    cpus 4

    input:
      set file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped), file(bestk) from rerun_admixture

    output:
      set file(vcf), file(vcfindex), val(nSM), val(maf), val(samples), val(ld), file(map), file(ped), file(bestk), file("log_BestK.out"), file("LD_${ld}_MAF_${maf}_BestK.P"), file("LD_${ld}_MAF_${maf}_BestK.Q") into admixture_bestk
      file("LD_${ld}_MAF_${maf}_BestK.Q") into admixture_bestk_to_popgenome

      """

        bk=`head -1 ${bestk}`

        admixture --cv=100 ${ped} \$bk -j4 | tee log_BestK.out

        mv LD_${ld}_MAF_${maf}.\$bk.P LD_${ld}_MAF_${maf}_BestK.P
        mv LD_${ld}_MAF_${maf}.\$bk.Q LD_${ld}_MAF_${maf}_BestK.Q
      """
}

/*
========================================
~ > *                              * < ~
~ ~ > *                          * < ~ ~
~ ~ ~ > *  HAPLOTYPE ANALYSIS  * < ~ ~ ~
~ ~ > *                          * < ~ ~
~ > *                              * < ~
========================================
*/

process_ibd=file("process_ibd.R")

/*
------------ Define IBDseq haplotype analysis parameters  
*/

minalleles = 0.02 // Specifies the minimum number of samples carrying the minor allele.
r2window = 1500 // Specifies the number of markers in the sliding window used to detect correlated markers.
ibdtrim = 0
r2max = 0.3

process ibdseq_haplotype {

    publishDir "${params.out}/HAPLOTYPE", mode: 'copy'

    cpus 8

    input:
      set file(vcf), file(vindex) from smallvcf_haplotype

    output:
      file("haplotype.tsv") into haplotype_analysis

      """

        minalleles=\$(bcftools query --list-samples ${vcf} | wc -l | awk '{ print \$0*${minalleles} }' | awk '{printf("%d\\n", \$0+=\$0<0?0:0.9)}')
        if [[ \${minalleles} -lt 2 ]];
        then
            minalleles=2;
        fi;
        echo "minalleles=${minalleles}"
        for chrom in I II III IV V X; do
            java -jar `which ibdseq.r1206.jar` \\
                gt=${vcf} \\
                out=haplotype_\${chrom} \\
                ibdtrim=${ibdtrim} \\
                minalleles=\${minalleles} \\
                r2max=${r2max} \\
                nthreads=${task.cpus} \\
                chrom=\${chrom}
            done;

        cat *.ibd | awk '{ print \$0 "\\t${minalleles}\\t${ibdtrim}\\t${r2window}\\t${r2max}" }' > haplotype.tsv
      """
}

/*
------------ Stitch together IBDseq segments  
*/

process analyze_ibdseq {

    publishDir "${params.out}/HAPLOTYPE", mode: 'copy'

    input:
        file("haplotype.tsv") from haplotype_analysis

    output:
        file("processed_haps.Rda")
        file("haplotype_plot_df.Rda") into plot_df


    """
      
      Rscript --vanilla `which process_ibd.R`
    """
}

/*
------------ Plot haplotypes, perform sweep analysis  
*/

process plot_ibdseq {

    publishDir "${params.out}/HAPLOTYPE", mode: 'copy'

    input:
        file("haplotype_plot_df.Rda") from plot_df

    output:
        file("haplotype_length.pdf")
        file("max_haplotype_sorted_genome_wide.pdf")
        file("haplotype.pdf")
        file("sweep_summary.tsv")

    """
      Rscript --vanilla `which plot_ibd.R`
    """
}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  Classical Population Genetics Statistics  * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

process popgenome_whole_pop {

  tag { CHROM }

  publishDir "${params.out}/POPGENOME/WINDOW/${CHROM}", mode: "copy"

  echo true

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_window
    each CHROM from contigs_popgenome_window

  output:
    set file("*Linkage_Statistics.Rda"), file("*Neutrality_Diversity_Statistics.Rda"), file("*WHOLE_POPULATION_Statistics.Rda") into popgenome_wholepop_statistics


  script:
    """
      #!/usr/bin/env Rscript

      require(PopGenome)
      require(data.table)
      require(tidyr)
      require(dplyr)
      require(glue)

      system(glue::glue("echo Initializing PopGenome Parameters"))

      chr1 <- c(1,15072434)
      chr2 <- c(1,15279421)
      chr3 <- c(1,13783801)
      chr4 <- c(1,17493829)
      chr5 <- c(1,20924180)
      chr6 <- c(1,17718942)
      
      chr.lengths <- list(chr1,chr2,chr3,chr4,chr5,chr6)
      chroms <- c("I","II","III","IV","V","X")
      
      ANALYSIS_CHROM <- "${CHROM}"
      CHROM_START <- chr.lengths[which(chroms == "${CHROM}")][[1]][1]
      CHROM_END <- chr.lengths[which(chroms == "${CHROM}")][[1]][2]

      SLIDE_DISTANCE <- ${params.popgenome_slide}
      WINDOW_SIZE <-  ${params.popgenome_window}

      OUTGROUP <- "${params.anc}"

      system(glue::glue("echo Done Initializing PopGenome Parameters - WindowSize = {WINDOW_SIZE}, StepSize = {SLIDE_DISTANCE}, Whole Population, Chromosome = ${CHROM}"))

      POPGENOME_VCF <- "${vcf}"
      POPGENOME_GFF <- "${params.gff}"

      system(glue::glue("echo PopGenome - Reading VCF file")) 
      GENOME_OBJECT <- PopGenome::readVCF(
                                          POPGENOME_VCF, 
                                          numcols = 10000, 
                                          tid = ANALYSIS_CHROM, 
                                          frompos = CHROM_START, 
                                          topos = CHROM_END, 
                                          approx = F)
     
      system(glue::glue("echo PopGenome - Setting Outgroup and Defining Window Size"))

      GENOME_OBJECT <- PopGenome::set.outgroup(GENOME_OBJECT, OUTGROUP,  diploid = FALSE)
      
      GENOME_OBJECT <- PopGenome::sliding.window.transform(
                                                            GENOME_OBJECT, 
                                                            width = WINDOW_SIZE, 
                                                            jump = SLIDE_DISTANCE,
                                                            type = 2, 
                                                            whole.data = FALSE
                                                          )

      system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Detail Stats"))
      GENOME_OBJECT <- PopGenome::detail.stats(GENOME_OBJECT)
      system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Neutrality Stats"))
      GENOME_OBJECT <- PopGenome::neutrality.stats(GENOME_OBJECT, detail = TRUE)
      system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Diversity Stats"))
      GENOME_OBJECT <- PopGenome::diversity.stats(GENOME_OBJECT, pi = TRUE)
      system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Fst Stats"))
      GENOME_OBJECT <- PopGenome::F_ST.stats(GENOME_OBJECT, mode = "nucleotide", detail = TRUE)
      system(glue::glue("echo PopGenome - Calculating Population Genetic Statistics - Linkage Stats"))
      GENOME_OBJECT <- PopGenome::linkage.stats(GENOME_OBJECT, do.ZnS = TRUE, do.WALL = TRUE)

      system(glue::glue("echo PopGenome - Finished Calculating Population Genetic Statistics - Saving File"))

      save(GENOME_OBJECT, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Statistics.Rda"))

      system(glue::glue("echo Generating PopGenome DataFrames - Extracting Window Bins"))

      windowStarts <- data.frame(snp_index = 1:length(colnames(GENOME_OBJECT@BIG.BIAL[[1]])),
                                 position = colnames(GENOME_OBJECT@BIG.BIAL[[1]]))
      
      slide_index <- cbind(data.frame(lapply(GENOME_OBJECT@SLIDE.POS,  function(x) as.numeric(floor(mean(x)))))) %>%
        tidyr::gather(temp, snp_index) %>%
        dplyr::select(-temp) %>%
        dplyr::left_join(., windowStarts, by = "snp_index")

      system(glue::glue("echo Generating PopGenome DataFrames - Extracting Linkage Stats"))

      linkage_df <- data.frame(Wall.B = c(GENOME_OBJECT@Wall.B),
                               Wall.Q = c(GENOME_OBJECT@Wall.Q),
                               Kelly.Z_nS = c(GENOME_OBJECT@Kelly.Z_nS),
                               Rozas.ZZ = c(GENOME_OBJECT@Rozas.ZZ),
                               Rozas.ZA = c(GENOME_OBJECT@Rozas.ZA)) %>%
                                 dplyr::mutate(Population = "WHOLE_POPULATION",
                                               WindowPosition = slide_index\$position) %>%
        tidyr::gather(LinkageStat, value, -Population, -WindowPosition)

      save(linkage_df, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Linkage_Statistics.Rda"))

      system(glue::glue("echo Generating PopGenome DataFrames - Extracting Neutrality and Diversity Stats"))

      neutrality_df <- data.frame(Tajima.D = c(GENOME_OBJECT@Tajima.D),
                                  n.segregating.sites = c(GENOME_OBJECT@n.segregating.sites),
                                  Fu.Li.F = c(GENOME_OBJECT@Fu.Li.F),
                                  Fu.Li.D = c(GENOME_OBJECT@Fu.Li.D),
                                  Fay.Wu.H = c(GENOME_OBJECT@Fay.Wu.H),
                                  Zeng.E = c(GENOME_OBJECT@Zeng.E),
                                  theta_Tajima = c(GENOME_OBJECT@theta_Tajima),
                                  theta_Watterson = c(GENOME_OBJECT@theta_Watterson),
                                  theta_Achaz.Watterson = c(GENOME_OBJECT@theta_Achaz.Watterson),
                                  theta_Achaz.Tajima = c(GENOME_OBJECT@theta_Achaz.Tajima),
                                  theta_Fay.Wu = c(GENOME_OBJECT@theta_Fay.Wu),
                                  theta_Zeng = c(GENOME_OBJECT@theta_Zeng),
                                  nuc.diversity.within = c(GENOME_OBJECT@nuc.diversity.within),
                                  PI = c(GENOME_OBJECT@Pi),
                                  theta_Fu.Li = c(GENOME_OBJECT@theta_Fu.Li),
                                  hap.diversity.within = c(GENOME_OBJECT@hap.diversity.within)) %>%
        dplyr::mutate(Population = "WHOLE_POPULATION",
                      WindowPosition = slide_index\$position) %>%
        tidyr::gather(statistic, value, -Population, -WindowPosition)
      
      save(neutrality_df, file = glue::glue("CHROMOSOME-{ANALYSIS_CHROM}_WHOLE_POPULATION_Neutrality_Diversity_Statistics.Rda"))
    """

}

process plot_popgenome_whole_pop {


  publishDir "${params.out}/POPGENOME/PLOTS", mode: "copy", pattern: "*.png"
  publishDir "${params.out}/POPGENOME/WHOLE_GENOME", mode: "copy", pattern: "*.Rda"

  input:
    file("*") from popgenome_wholepop_statistics.collect()

  output:
    set file("Ce_Genome-wide_Neutrality_stats.Rda"), file("Ce_Genome-wide_Linkage_stats.Rda") into whole_genome_popgenome
    file("*.png") into whole_genome_popgenome_plots


  script:

  """
    #!/usr/bin/env Rscript
    library(tidyverse)

    # PLOT NEUTRALITY STATISTICS
    Neutrality_files <- grep(glue::glue("Neutrality"),list.files(), value = T)
    
    neutrality_ls <- list()
    for( i in 1:length(Neutrality_files) ){
      
      chr <- strsplit(strsplit(Neutrality_files[i],split = "_")[[1]][1], split="-")[[1]][[2]]
      if(chr == "23"){
        chr <- "X"
      }

      load(Neutrality_files[i])
      neutrality_ls[[i]] <- data.frame(neutrality_df) %>%
        dplyr::mutate(CHROM = chr)
    }

    neutrality_df <- dplyr::bind_rows(neutrality_ls)
    
    save(neutrality_df, file = glue::glue("Ce_Genome-wide_Neutrality_stats.Rda"))

    for(neutrality_statistic in unique(neutrality_df\$statistic)){

      neutrpl <- neutrality_df %>%
        dplyr::filter(statistic == neutrality_statistic) %>%
        ggplot()+
        aes(x = WindowPosition/1e6, y = value) +
        geom_point(alpha = 0.25 ) +
        geom_point(alpha = 0.5, size = 0.5) +
        theme_bw()+
        facet_grid(.~CHROM, scales ="free", space = "free_x") +
        theme(axis.text.x = ggplot2::element_text(size = 16),
              axis.text.y = ggplot2::element_text(size = 16),
              axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
              axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3)) +
        labs(x = "Genomic Position (Mb)",
             y = neutrality_statistic)
      
      ggsave(plot = neutrpl,
             glue::glue("Genome-wide_{neutrality_statistic}.png"),
             width = 12,
             height = 4,
             dpi = 300)
    }


    Linkage_files <- grep(glue::glue("Linkage"),list.files(), value = T)
  
    linkage_ls <- list()
    for( i in 1:length(Linkage_files)){
      
      chr <- strsplit(strsplit(Linkage_files[i],split = "_")[[1]][1], split="-")[[1]][[2]]
      if(chr == "23"){
        chr <- "X"
      }
      load(Linkage_files[i])
      linkage_ls[[i]] <- data.frame(linkage_df) %>%
        dplyr::mutate(CHROM = chr)
    }
    
    linkage_df <- dplyr::bind_rows(linkage_ls)
    
    save(linkage_df, file = glue::glue("Ce_Genome-wide_Linkage_stats.Rda"))
    
    for(linkage_statistic in unique(linkage_df\$LinkageStat)){

      linkpl <- linkage_df %>%
        dplyr::filter(LinkageStat == linkage_statistic) %>%
        ggplot()+
        aes(x = WindowPosition/1e6, y = value, color = Population) +
        geom_point(size = 0.5, alpha = 0.5)+
        theme_bw() +
        facet_grid(.~CHROM, scales ="free", space = "free_x") +
        theme(axis.text.x = ggplot2::element_text(size = 16),
              axis.text.y = ggplot2::element_text(size = 16),
              # legend.position = "none",
              axis.title.x = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3),
              axis.title.y = ggplot2::element_text(size = 18, face = "bold", color = "black", vjust = -0.3)) +
        labs(x = "Genomic Position (Mb)",
             y = linkage_statistic)
      
      ggsave(plot = linkpl, 
             glue::glue("Genome-wide_{linkage_statistic}.png"),
             width = 12,
             height = 4,
             dpi = 300)
    }

  """


}

/*
==============================================================
~ > *                                                    * < ~
~ ~ > *                                                * < ~ ~
~ ~ ~ > *  Classical PopGen Statistics - Gene Level * < ~ ~ ~
~ ~ > *                                                * < ~ ~
~ > *                                                    * < ~
==============================================================
*/

/*
------------ Calculate Gene-level PopGen Statistics
*/

process popgenome_gene_complete {

  tag { "${CHROM}" }

  publishDir "${params.out}/POPGENOME/GENE/${CHROM}", mode: "copy"

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_gene
    each CHROM from contigs_popgenome_gene

  output:
    set val("${CHROM}"), file("CHROMOSOME-${CHROM}_WHOLE_POPULATION_Statistics.Rda"), file("CHROMOSOME-${CHROM}_WHOLE_POPULATION_PopGenome_Object.Rda") into fst_complete_gene_output

    """

      Rscript --vanilla `which PopGenome_Gene_exons.R` ${CHROM} ${params.anc} ${vcf} ${params.gff}
    """

}

process popgenome_gene_kpop {

  tag { "${CHROM}" }

  publishDir "${params.out}/POPGENOME/SUBPOPs/GENE/${CHROM}", mode: "copy"

  input:
    set file(vcf), file(vindex) from smallvcf_classic_popgen_gene_admix_pops
    file(popfile) from admixture_bestk_to_popgenome
    each CHROM from contigs_popgenome_gene_subpops

  output:
    set val("${CHROM}"), file("CHROMOSOME-${CHROM}_WHOLE_POPULATION_FST_Statistics.Rda"), file("CHROMOSOME-${CHROM}_WHOLE_POPULATION_ND_Statistics.Rda") into fst_pops_gene_output

    """

      Rscript --vanilla `which PopGenome_Gene_exons_subpops.R` ${CHROM} ${params.anc} ${vcf} ${params.gff} ${popfile}
    """

}

/*
======================================================================
~ > *                                                            * < ~
~ ~ > *                                                        * < ~ ~
~ ~ ~ > *  Fancier Selection Population Genetics Statistics  * < ~ ~ ~
~ ~ > *                                                        * < ~ ~
~ > *                                                            * < ~
======================================================================
*/



process vcf_to_h21 {

  publishDir "${params.out}/H21/INPUT/", mode: "copy"

  tag {CHROM}

  input:
    set file(vcf), file(vindex) from smallvcf_h21
    each CHROM from contigs_h21

  output:
    set val(CHROM), file("${CHROM}_H21_input.txt") into h21_input

  """

    bcftools view -v snps ${vcf} |\\
    bcftools query --regions ${CHROM} -f '%POS,[%TGT,]\\n' |\\
    sed 's/,\$//' |\\
    sed -E 's/(.)\\/(.)/\\2/g' |\\
    sed 's/\\./N/g' > ${CHROM}_H21_input.txt
  """

}

// this runs H21 from https://github.com/ngarud/SelectionHapStats
// still need simulations via msms to rule out neutrality

/*process run_h21 {

  publishDir "${params.out}/H21/OUTPUT/", mode: "copy"

  tag {CHROM}

  input:
    set val(CHROM), file(h21_in) from h21_input

  output:
    file("${CHROM}_H21_Results.txt") into h21_output

  """

    samples=`head -1 ${h21_in} | tr -cd ',' | wc -c`
    python2 `which H12_H2H1.py` ${h21_in} \$samples -o ${CHROM}_H21_Results.txt -w 200 -j 25 -d 0
  """
}*\


Channel.from(0.003, 0.006, 0.012, 0.024, 0.048, 0.1)
       .spread(smallvcf_raisd)
       .into{raisd_input}


process run_raisd {

  tag {"${CHROM} - ${maf}"}

  publishDir "${params.out}/RAISD/${maf}", mode: "copy"

  input:
    set val(maf), file(vcf), file(vindex) from raisd_input
    each CHROM from contigs_raisd

  output:
    set val("${maf}"), file("RAiSD_Info*"), file("RAiSD_Report*") into raisd_output

    """

    bcftools view --regions=${CHROM} -Ov -o ${CHROM}_raisd_input.vcf ${vcf}

    ${params.raisd} \\
    -n ${CHROM}_${maf} \\
    -I ${CHROM}_raisd_input.vcf \\
    -m ${maf}
    """
}

/*
======================================
~ > *                            * < ~
~ ~ > *                        * < ~ ~
~ ~ ~ > *  Run PCA and DAPC  * < ~ ~ ~
~ ~ > *                        * < ~ ~
~ > *                            * < ~
======================================
*/

/*
------------ Prepare files for EIGENSTRAT
*/

process vcf_to_eigstrat_files {

    tag {"PREPARE EIGENSTRAT FILES"}

    publishDir "${params.out}/EIGESTRAT/INPUTFILES/", mode: 'copy'

    input:
      set file(vcf), file(vcfindex) from annotated_vcf_eigenstrat

    output:
      set file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind") into eigenstrat_input

      """

      bcftools view --regions I,II,III,IV,V,X ${vcf} |\\
      bcftools norm -m +snps -Oz -o ce_norm.vcf.gz

      tabix -p vcf ce_norm.vcf.gz

      plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 0.95 --allow-extra-chr 

      plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recode12 --out eigenstrat_input --allow-extra-chr

      awk -F":" '\$1=\$1' OFS="\\t" plink.prune.in | \\
      sort -k1,1d -k2,2n > markers.txt

      bcftools query -l ce_norm.vcf.gz |\\
      sort > sorted_samples.txt 

      bcftools view -v snps -S sorted_samples.txt -R markers.txt ce_norm.vcf.gz |\\
      bcftools query -f '%CHROM\\t%CHROM:%POS\\t%cM\\t%POS\\t%REF\\t%ALT\\n' ce_norm.vcf.gz |\\
      sed 's/^III/3/g' |\\
      sed 's/^II/2/g' |\\
      sed 's/^IV/4/g' |\\
      sed 's/^I/1/g' |\\
      sed 's/^V/5/g' > eigenstrat_input.pedsnp      

      cut -f-6 -d' ' eigenstrat_input.ped |\\
      awk '{print 1, \$2, \$3, \$3, \$5, 1}'  > eigenstrat_input.pedind

      echo "rerun"
      """

}

eigenstrat_input
  .into{eigenstrat_no_outlier;
        eigenstrat_outlier_removal;
        eigenstrat_fst}


/*
------------ Run EIGENSTRAT without removing outlier strains
*/

process run_eigenstrat_no_outlier_removal {

    publishDir "${params.out}/EIGESTRAT/NO_REMOVAL/", mode: 'copy'

    input:
      set file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind") from eigenstrat_no_outlier
      file(eigenparameters) from eigenstrat_noremoval

    output:
      set file("eigenstrat_no_removal.evac"), file("eigenstrat_no_removal.eval"), file("logfile_no_removal.txt"), file("eigenstrat_no_removal_relatedness"), file("eigenstrat_no_removal_relatedness.id"), file("TracyWidom_statistics_no_removal.tsv") into eigenstrat_outlier_removal_output

      """

      smartpca -p ${eigenparameters} > logfile_no_removal.txt

      sed -n -e '/Tracy/,\$p' logfile_no_removal.txt |\
      sed -e '/kurt/,\$d' |\
      awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
      sed -e "s/[[:space:]]\\+/ /g" |\
      sed 's/^ //g' |\
      awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
      awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_no_removal.tsv
      """

}

/*
------------ Run EIGENSTRAT with removing outlier strains
*/

process run_eigenstrat_with_outlier_removal {


    publishDir "${params.out}/EIGESTRAT/OUTLIER_REMOVAL/", mode: 'copy'

    input:
      set file("eigenstrat_input.ped"), file("eigenstrat_input.pedsnp"), file("eigenstrat_input.pedind") from eigenstrat_outlier_removal
      file(eigenparameters) from eigenstrat_removal

    output:
      set file("eigenstrat_outliers_removed.evac"), file("eigenstrat_outliers_removed.eval"), file("logfile_outlier.txt"), file("eigenstrat_outliers_removed_relatedness"), file("eigenstrat_outliers_removed_relatedness.id"), file("TracyWidom_statistics_outlier_removal.tsv") into eigenstrat_no_outlier_removal_output
      
      """

      smartpca -p ${eigenparameters} > logfile_outlier.txt

      sed -n -e '/Tracy/,\$p' logfile_outlier.txt |\
      sed -e '/kurt/,\$d' |\
      awk '\$0 !~ "##" && \$0 !~ "#" {print}' |\
      sed -e "s/[[:space:]]\\+/ /g" |\
      sed 's/^ //g' |\
      awk 'BEGIN{print "N", "eigenvalue", "difference", "twstat", "p-value", "effect.n"}; {print}' OFS="\\t" |\
      awk -F" " '\$1=\$1' OFS="\\t" > TracyWidom_statistics_outlier_removal.tsv
      """

}

/*
------------ Prepare files for DAPC analysis
*/

process prep_files_for_DAPC {


    publishDir "${params.out}/DAPC/INPUT_FILES", mode: 'copy'

    input:
      set file(vcf), file(vcfindex) from smallvcf_dapc

    output:
      file("dapc_input.raw") into dapc_input
      
      """

      bcftools view --regions I,II,III,IV,V,X ${vcf} |\\
      bcftools norm -m +snps -Oz -o ce_norm.vcf.gz

      tabix -p vcf ce_norm.vcf.gz

      plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --indep-pairwise 50 10 0.95 --allow-extra-chr 

      plink --vcf ce_norm.vcf.gz --snps-only --biallelic-only --set-missing-var-ids @:# --extract plink.prune.in --geno --recodeA --out dapc_input --allow-extra-chr
      """

}


Channel.from( 3,4,5,6,7,8 )
       .spread(dapc_input)
       .into{dapc_pop_input}

/*
process run_dapc {

  tag {"${pops}"}

  cpus 8

  publishDir "${params.out}/DAPC/RESULTS/${pops}", mode: "copy"

  input:
    set val(pops), file(genotypes) from dapc_pop_input

  output:
    set val("DAPC_LOADING.pdf"), file("DAPC_SCATTER.pdf"), file("Adegent_DAPC_Object.Rda"), file("Adegent_PC_Object.Rda"), file("SNP_LOADINGS.tsv"), file("STRAIN_MEMBERSHIP_PROBABILITIES.tsv") into dapc_output

    """
      Rscript --vanilla `which run_dapc.R` ${genotypes} ${pops} ${task.cpus}
    """

}
*/