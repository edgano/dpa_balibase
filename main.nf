#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2019, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'regressive-msa-analysis'.
 *
 *   regressive-msa--analysis is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   regressive-msa-analysis is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with regressive-msa-analysis.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main regressive-msa-analysis pipeline script
 *
 * @authors
 * Evan Floden <evanfloden@gmail.com>
 * Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * Edgar Garriga
 * Cedric Notredame 
 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
params.tfa = "$baseDir/data/RV11/*.tfa"

// input reference sequences aligned in 
params.refs = "$baseDir/data/RV11/*.msf"

// input guide trees in Newick format. Or `false` to generate trees
//params.trees = "$baseDir/data/trees/*.CLUSTALO.dnd"
params.trees = false

// which alignment methods to run
params.align_method = "CLUSTALO,MAFFT-FFTNS1" //"UPP,MAFFT-SPARSECORE,MAFFT-FFTNS1,CLUSTALO,MAFFT-GINSI,PROBCONS,UPP"

// which tree methods to run if `trees` == `false`
params.tree_method = "CLUSTALO" //CLUSTALO-RANDOM,CLUSTALO,MAFFT-FFTNS1,MAFFT_PARTTREE"

// generate regressive alignments ?
params.regressive_align = true

// create standard alignments ?
params.standard_align = false

// create default alignments ? 
params.default_align = false

// evaluate alignments ?
params.evaluate = true

// bucket sizes for regressive algorithm
params.buckets= '10,50,100,500,1000'

// output directory
params.output = "$baseDir/results"

params.reformat = true

log.info """\
         R E G R E S S I V E   M S A   A n a l y s i s  ~  version 0.1"
         ======================================="
         Input TFS       (TFA)                          : ${params.tfa}
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA)               : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Output directory (DIRECTORY)                   : ${params.output}
         Alignment methods                              : ${params.align_method}
         Tree methods                                   : ${params.tree_method}
         Reformat TFA to FA                             : ${params.reformat}
         Generate default alignments                    : ${params.default_align}
         Generate standard alignments                   : ${params.standard_align}
         Generate regressive alignments (DPA)           : ${params.regressive_align}
         Bucket Sizes for regressive alignments         : ${params.buckets}
         Perform evaluation? Requires reference         : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.output}
         \
         """
         .stripIndent()


// Channels containing sequences
if ( params.tfa ) {
  Channel
  .fromPath(params.tfa)
  .map { item -> [ item.baseName, item] }
  .set { tfaSeqs }
}

// Channels containing reference alignments for evaluation [OPTIONAL]
if( params.refs ) {
  Channel
    .fromPath(params.refs)
    .map { item -> [ item.baseName, item] }
    .set { refs }
}

// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
    .set { trees }
}
else { 
  Channel
    .empty()
    .set { trees }
}

tree_methods = params.tree_method
align_methods = params.align_method

process reformatBaliBase {
    tag "${id} - ${seqs}"
    publishDir "${params.output}/data", mode: 'copy', overwrite: true

    input:
      set val(id), \
         file(seqs) \
         from tfaSeqs

    output:
     set val(id), \
       file("${id}.fa") \
       into reformatSeqs, reformatSeqs2

    when:
      params.reformat

    script:
      """
      t_coffee -other_pg seq_reformat -in ${seqs} -output=fasta_seq -out=${id}.fa
      """
}
process guide_trees {

    tag "${id}.${tree_method}"
    publishDir "${params.output}/guide_trees", mode: 'copy', overwrite: true
   
    input:

     set val(id), \
         file(seqs) \
         from reformatSeqs
     each tree_method from tree_methods.tokenize(',') 

   output:
     set val(id), \
       val(tree_method), \
       file("${id}.${tree_method}.dnd") \
       into treesGenerated

   when:
     !params.trees

   script:
     """
     clustalo -i ${seqs} --guidetree-out ${id}.${tree_method}.dnd --force
     """
}

treesGenerated
  .mix ( trees )
  .combine ( reformatSeqs2, by:0 )
  .into { seqsAndTreesForStandardAlignment; seqsAndTreesForRegressiveAlignment }

process standard_alignment {
  
    tag "${id}.${align_method}.STD.NA.${tree_method}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true
    
    input:
      set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForStandardAlignment

      each align_method from align_methods.tokenize(',') 

    when:
      params.standard_align

    output:
      set val(id), \
      val("${align_method}"), \
      val(tree_method), val("std_align"), \
      val("NA"), \
      file("*.aln"), \
      file ("*.msf") \
      into standard_alignments

     script:
     """
     clustalo --infile=${seqs} \
        --guidetree-in=${guide_tree} \
        --outfmt=fa \
        -o ${id}.std.${align_method}.with.${tree_method}.tree.aln
      
      t_coffee -other_pg seq_reformat -in ${id}.std.${align_method}.with.${tree_method}.tree.aln -output=msf_aln -out=${id}.std.${align_method}.with.${tree_method}.tree.aln.msf
     """
}

process regressive_alignment {

    tag "${id}.${align_method}.DPA.${bucket_size}.${tree_method}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input:
      set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForRegressiveAlignment

      each bucket_size from params.buckets.tokenize(',')
       
      each align_method from align_methods.tokenize(',')   

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("dpa_align"), \
        val(bucket_size), \
        file("*.aln"), \
        file ("*.msf") \
        into regressive_alignments , regressive_alignments2

    when:
      params.regressive_align

    script:
    """
    t_coffee -reg -reg_method clustalo_msa \
         -reg_tree ${guide_tree} \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -outfile ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln
    
    t_coffee -other_pg seq_reformat -in ${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tree.aln -output=msf_aln -out=${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.aln.msf
    """
}

standard_alignments
  .mix ( regressive_alignments )
  .set { all_alignments }

refs
  .cross (all_alignments)
  .map { it -> [it[0][0], it[1][1], it[1][2], it[1][4], it[1][6], it[0][1]] }
  .set { toEvaluate }

process eval{
  
  tag "${id}"
  publishDir "${params.output}/individual_scores", mode: 'copy', overwrite: true
    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(bucket_size), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set val(id), \
          val(tree_method), \
          val(align_method), \
          val(bucket_size), \
          file("*.sp"), \
          file("*.tc") \
          into scores

    when:
      params.evaluate

    script:
    """
    ${baseDir}/bin/bali_score ${ref_alignment} ${test_alignment} -v > score.out
    awk -F " " 'END{print \$3 > "${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.sp"; print \$4 > "${id}.dpa_${bucket_size}.${align_method}.with.${tree_method}.tc" }' score.out 
    """

}
workflow.onComplete {
  // println (['bash','-c', "$baseDir/bin/cpu_calculate.sh ${params.output}/individual_scores"].execute().text)
  println "Execution status: ${ workflow.success ? 'OK' : 'failed' } runName: ${workflow.runName}" 
}