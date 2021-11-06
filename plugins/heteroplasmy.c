/*  plugins/heteroplasmy.c -- Heteroplamic Mitochondrial Variant Caller

    Copyright (C) 2021 Grant Daly

    Author: Grant Daly <daly@southalabama.edu>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>

#include <stdlib.h>

#include <inttypes.h>

#include <htslib/vcf.h>

static bcf_hdr_t * hdr;

/*
    This short description is used to generate the output of `bcftools plugin -l`.
*/
const char * about(void) {
  return "Heteroplasmic Mitochondrial Variant Caller\n";

}

/*
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/
int init(int argc, char ** argv, bcf_hdr_t * in , bcf_hdr_t * out) {
  const char vafHead[] = "##FORMAT=<ID=AF,Number=R,Type=Float,Description=\"Allele Frequency\"";

  //out = in;
  // if the AF tag isn't present, add it
  if (bcf_hdr_id2int(out, BCF_DT_ID, "AF") == -1) {
    if (bcf_hdr_append(out, vafHead) == 0) {
      //printf("appended\n");
    }  else {
      printf("failed append\n");
      exit(1);
    }

  }

  const char artHead[] = "##FORMAT=<ID=ART,Number=A,Type=Integer,Description=\"Artifact Score\"";
  // if the AF tag isn't present, add it
  if (bcf_hdr_id2int(out, BCF_DT_ID, "ART") == -1) {
    if (bcf_hdr_append(out, artHead) == 0) {
      //printf("appended\n");
    } else {
      printf("failed append\n");
      exit(1);
    }

  }
  // hdr is the static header I will be accessing
  hdr = out;

  return 0;
}

typedef struct sample_param {
  int n_allele;
  int dp;
int sampleIndex;
float *afArray;
int * adArray;
int * gtArray;
  int * adfArray;
int * adrArray;
  int * artArray;
  }
sample_param_t;

int minimumInt(int a, int b){
  if(a <= b)
    return a;
  else
    return b;
}

int maximumInt(int a, int b){
  if(a >= b)
    return a;
  else
    return b;
}
void processSample(sample_param_t params) {
  
  const int n_allele = params.n_allele;
  //const int tempDP = params.dp;

  
  float * afArray = params.afArray;
  int * adArray = params.adArray;
  int * adfArray = params.adfArray;
  int * adrArray = params.adrArray;

  // check if sample with missing values. Using first val of adArray
  
  if(adArray[0] == bcf_int32_missing) {
    return;
  }
  else if(adArray[0] < 0) {
    return;
  }
  int tempAD = 0;
  int tempADF = 0;
  int tempADR = 0;
  float tempAF = 0.0;
  
  int * artArray = params.artArray;

  int sampleIndex = params.sampleIndex;

  int * gt_arr = params.gtArray;



  // call possible artifacts on non reference alleles
  int tempADFRef = adfArray[(sampleIndex * n_allele)];
  int tempADRRef = adrArray[(sampleIndex * n_allele)];
  int tempADRef = tempADFRef + tempADRRef;
  int tempADFVar = 0;
  int tempADRVar = 0;
  int numAltAlleles = n_allele - 1;

  if((tempADRef >= 100)){
      int minimumADStrand =  minimumInt(tempADFRef, tempADRRef);
      if(tempADRef != 0){
	double refFrac = minimumADStrand / ((double) tempADRef);

	// if meets the minimums so calculate artifact score
	if(refFrac >= 0.2){
	  for (int j = 1; j < n_allele ; j++){
	    tempADFVar = adfArray[(sampleIndex * n_allele) + j];
	    tempADRVar = adrArray[(sampleIndex * n_allele) + j];
	    //printf("forward %d reverse %d\n", tempADFVar, tempADRVar);
	    if((tempADFVar >= 3) && (tempADRVar == 0)) {
	      //printf("forward artifact\n");
	      artArray[(sampleIndex * numAltAlleles) + (j - 1)] = tempADFVar;
	      }
	    else if((tempADRVar >= 3) && (tempADFVar == 0)) {
	      //printf("reverse artifact\n");
	      artArray[(sampleIndex * numAltAlleles) + (j - 1)] = - 1 * tempADRVar;
	      }
	    //printf("array value %d\n", artArray[(sampleIndex * numAltAlleles) + (j - 1)]);
	  }
      }
  }
  
    

  }
  
  
  int candGT_arr[2] = {
      -1,
      -1
    };
  //ndGT_arr[0] = -1;
  //candGT_arr[1] = -1;

  int totalBases = 0;
  for (int j = 0; j < n_allele; j++) {
    totalBases += adArray[(sampleIndex * n_allele) + j];
  }
  // call variants on all alleles

  // normal vcf has 2 alleles so I'm using this. pretty good assumption because there are unlikely to be many tripplet heteroplasmies
  //int maxAlt = 2;
  //int numberAllelesAssigned = 0;

  //calculate the allele frequencies. I am separating this because the variant calling logic has some breaks when a condition is met. This meant the AF calculation was sometimes scipped
  for (int j = 0; j < n_allele; j++) {
    tempAD = adArray[(sampleIndex * n_allele) + j];
    tempADF = adfArray[(sampleIndex * n_allele) + j];
    tempADR = adrArray[(sampleIndex * n_allele) + j];

    // calculate AF
    
    if (totalBases > 0) {

      tempAF = (float) tempAD / (float) totalBases;
      //printf("allele depth %d total bases %d fraction %f\n",tempAD, totalBases, tempAF);
    } else {
      tempAF = 0;
    }
    afArray[(sampleIndex * n_allele) + j] = tempAF;

  }


  
  for (int j = 0; j < n_allele; j++) {
    tempAD = adArray[(sampleIndex * n_allele) + j];
    tempADF = adfArray[(sampleIndex * n_allele) + j];
    tempADR = adrArray[(sampleIndex * n_allele) + j];

    if((totalBases >=20) && (tempADF >= 10) && (tempADR >= 10)) {

      // retreiving the AF for this allele from the array we previously filled
    tempAF = afArray[(sampleIndex * n_allele) + j];
    if (tempAF >= 0.99) {
      candGT_arr[0] = j;
      candGT_arr[1] = j;
      break;
      // if I find a homozygous site can skip additional searching
      // could put a goto if I wanted minor improvement
    }  else if ((totalBases >= 500) && (tempAF >= 0.01)) {
      // heteroplasmy call
      if (candGT_arr[0] == -1) {

        candGT_arr[0] = j;
      }
      else if (candGT_arr[1] == -1) {
        candGT_arr[1] = j;
      }
      else {
	break;
        //	  printf("site has more that 2 heteroplamies, and cannot display additional\n");
      }
    }
  }
  }
  //printf("slot 1 %d slot 2 %d\n", candGT_arr[0], candGT_arr[1]);
    // verify that these values were changed, and if so update genotype
    if (candGT_arr[0] != -1){
      // doing order in reverse to prefer 0/1 opposed to 1/0.
      gt_arr[sampleIndex * 2 + 1] = bcf_gt_unphased(candGT_arr[0]);
    }
    else
      {
	gt_arr[sampleIndex * 2 + 1] = bcf_gt_unphased(0);
      }
    if (candGT_arr[1] != -1){
      gt_arr[sampleIndex * 2] = bcf_gt_unphased(candGT_arr[1]);
    }
    else
      {
	gt_arr[sampleIndex * 2] = bcf_gt_unphased(0);
      }
    

  

}
/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t * process(bcf1_t * rec) {

  bcf_unpack(rec, BCF_UN_ALL);

  const int n_allele = rec -> n_allele;
  const int nsamples = bcf_hdr_nsamples(hdr);

  int32_t * adArray = 0;
  int adCount = 0;
  //int formatReturn;
  //int formatReturn = bcf_get_format_int32(hdr, rec, "AD", adArray, adCount);
  if (bcf_get_format_int32(hdr, rec, "AD", & adArray, & adCount) < 0) {
    printf("error with get format\n");
    exit(1);
  }

  int32_t * adfArray = 0;
  int adfCount = 0;
  if (bcf_get_format_int32(hdr, rec, "ADF", & adfArray, & adfCount) < 0) {
    printf("error with get format\n");
    exit(1);
    }

  int32_t * adrArray = 0;
  int adrCount = 0;
  if (bcf_get_format_int32(hdr, rec, "ADR", & adrArray, & adrCount) < 0) {
    printf("error with get format\n");
    exit(1);
    }

  // get dp (number of high quality bases)
  /* int32_t * depthDP = 0; */
  /* int dpCount = 0; */
  /* if (bcf_get_format_int32(hdr, rec, "DP", & depthDP, & dpCount) < 0) { */
  /*   printf("error with get format\n"); */
  /*   exit(1); */
  /* } */

  // get genotypes
  int ngt, * gt_arr = NULL, ngt_arr = 0;
  ngt = bcf_get_genotypes(hdr, rec, & gt_arr, & ngt_arr);
  //printf("ngt %d    ngt_arr %d\n", ngt, ngt_arr);
  /*for(int i=0; i < ngt;i++){
    printf("genotypes %d\n", gt_arr[i]);
    }*/

  // get allele frequency array
  /*double *afArray = NULL;
  int afCount = 0;
  int afReturn = bcf_get_format_float(hdr, inLocus, "DP", &afArray, &afCount);
  printf("######## number of af %d or %d\n", afReturn, afCount);*/
  int numAF = n_allele * nsamples;
  float * afArray = malloc(numAF * sizeof(float));
  for(int i = 0; i < numAF; i++){
    //afArray[i] = bcf_float_missing;
    bcf_float_set_missing(afArray[i]);
  }


  // create and initialize the artifact array.
  int numART = (n_allele -1) * nsamples;
  int * artArray = malloc(numART * sizeof(int));
  for(int i = 0; i < numART; i++){
    artArray[i] = bcf_int32_missing;
  }

  //int32_t tempAD = 0;
  //int32_t tempDP = 0;
  for (int i = 0; i < nsamples; i++) {
    // check if the genotype is missing from low coverage sites and skip if soft
    // causes issue because AF ends up not being initialized. could possibly put this after the AF logic.
    /*if (bcf_gt_is_missing(gt_arr[i*2])) {
	//	printf("############ missing gt ############\n");
	continue;
	}*/

    /* int candGT_arr[2] = { */
    /*   0, */
    /*   0 */
    /* }; */

    //tempDP = depthDP[i];

    
    sample_param_t sampleParams = {
      .n_allele = n_allele,
      .gtArray = gt_arr,
      .adArray = adArray,
      .afArray = afArray,
      .sampleIndex = i,
      .adfArray = adfArray,
      .adrArray = adrArray,
      .artArray = artArray
      
    };
    //processSample(n_allele, tempDF, tempAD, afArray);
    processSample(sampleParams);
    
  }

  bcf_update_genotypes(hdr, rec, gt_arr, ngt);

  // update format AF tags
  bcf_update_format_float(hdr, rec, "AF", afArray, numAF);
  
  if(bcf_update_format_int32(hdr, rec, "ART", artArray, numART) != 0){
    printf("error updating artifact score\n");
  }

  free(afArray);
  free(artArray);
  return rec;

}

/*
    Clean up.
*/
void destroy(void) {

}
