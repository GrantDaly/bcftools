/*  plugins/heteroplasmy.c -- Heuristic Mitochondrial Variant Caller

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

static bcf_hdr_t *hdr;


/*
    This short description is used to generate the output of `bcftools plugin -l`.
*/
const char *about(void)
{
    return
      "Heuristic Mitochondrial Variant Caller\n";
        
}


/*
    Called once at startup, allows to initialize local variables.
    Return 1 to suppress VCF/BCF header from printing, 0 otherwise.
*/    
int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{
    char vafHead[] = "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\"";

    //out = in;
  // if the AF tag isn't present, add it
  if(bcf_hdr_id2int(out, BCF_DT_ID, "AF") == -1) {
    if( bcf_hdr_append(out, vafHead) == 0){ 
      //printf("appended\n");
    }
    else
    {
      printf("failed append\n");
      exit(1);
    }
    
  }
  // hdr is the static header I will be accessing
  hdr = out;
    
    return 0;
}


/*
    Called for each VCF record. Return rec to output the line or NULL
    to suppress output.
*/
bcf1_t *process(bcf1_t *rec)
{

      bcf_unpack(rec, BCF_UN_ALL);

      const int n_allele = rec->n_allele;
      const int nsamples = bcf_hdr_nsamples(hdr);
    
    int32_t *values = 0;
    int count = 0;
    int formatReturn = bcf_get_format_int32(hdr, rec, "AD", &values, &count);
    if(formatReturn < 0){
      printf("error with get format\n");
      exit(1);
    }

    // get dp (number of high quality bases)
    int32_t *depthDP = 0;
    int dpCount = 0;
    formatReturn = bcf_get_format_int32(hdr, rec, "DP", &depthDP, &dpCount);
    
    // get genotypes
    int ngt, *gt_arr = NULL, ngt_arr = 0;
    ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
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
    float* afArray = malloc(n_allele * nsamples * sizeof(float));
    
    int32_t tempAD = 0;
    int32_t tempDP = 0;
    for(int i =0; i < nsamples;i++){
      // check if the genotype is missing from low coverage sites and skip if soft
      // causes issue because AF ends up not being initialized. could possibly put this after the AF logic.
      /*if (bcf_gt_is_missing(gt_arr[i*2])) {
	//	printf("############ missing gt ############\n");
	continue;
	}*/
      
      int candGT_arr[2] = {0, 0};
      
      tempDP = depthDP[i];
      
      candGT_arr[0] = -1;
      candGT_arr[1] = -1;
      for(int j =0; j < n_allele; j++){
      
      tempAD = values[(i*n_allele) + j];
      
      float tempAF = 0.0;
      if(tempDP != 0){
	
	tempAF = (float)((double) tempAD / (double) tempDP) ;
      }
      else
	{
	  tempAF = 0;
	}
      
      afArray[(i*n_allele) + j] = tempAF;

      if((tempDP >= 20) && (tempAF >= 0.99)){
	candGT_arr[0] = j;
	candGT_arr[1] = j;
	// if I find a homozygous site can skip additional searching
	// could put a goto if I wanted minor improvement
	}
      else if(tempDP < 20){
	// if depth requirment not met make it the original values. can be missing genotype
      }
      
      else if((tempDP >= 500) && (tempAF >= 0.01)){
	// heteroplasmy call, need to see if the pigeon hole is filled
	if(candGT_arr[0] == -1){
 	  candGT_arr[0] = j;
	}
	else if(candGT_arr[1] == -1){
	  candGT_arr[1] = j;
	}
	else {
	  //	  printf("site has more that 2 heteroplamies, and cannot display additional\n");
      }
      }

      }
      
      // verify that these values were changed, and if so update genotype
      if(candGT_arr[0] != -1)
      gt_arr[i*2] = bcf_gt_unphased(candGT_arr[0]);

      if(candGT_arr[1] != -1)
      gt_arr[i*2 + 1] = bcf_gt_unphased(candGT_arr[1]);

    }

    bcf_update_genotypes(hdr, rec, gt_arr, ngt);

    // update format AF tags
    bcf_update_format_float(hdr, rec, "AF", afArray, numAF);
    
     free(afArray);
     return rec;

}


    


/*
    Clean up.
*/
void destroy(void)
{

}


