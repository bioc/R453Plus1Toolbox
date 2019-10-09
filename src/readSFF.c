#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "io_lib/os.h"
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include "readSFF.h"
#include <R.h>
#include <Rdefines.h>

/***************************************************
 * Function to read in a Roche .sff file creating  *
 * the container structure described in the header *
 * consisting of one header and many reads.        *
 ***************************************************/

sff_container *readSFF(char *filename) {

  /**********************************************
   * Read the given file completely into buffer *
   **********************************************/

  FILE *file;
  char *buffer;
  unsigned long fileLen;
  int block_count;

  if (NULL == (file = fopen(filename, "rb"))) {
    REprintf("Unable to open file %s \n", filename);
    return NULL;
  }

  //Get file length
  fseek(file, 0, SEEK_END);
  fileLen = ftell(file);
  fseek(file, 0, SEEK_SET);

  //Allocate memory
  if (NULL == (buffer = (char *)malloc(fileLen+1))) {
    REprintf("Memory error!\n");
    fclose(file);
    return NULL;
  }

  //Read file contents into buffer
  block_count = fread(buffer, fileLen, 1, file);
  fclose(file);


  /******************************
   * Read in the header section *
   ******************************/

  //Pointer for position in buffer
  int pointer;

  sff_header *header;
  if(NULL == (header = (sff_header *)calloc(1, sizeof(*header)))) {
    REprintf("Memory error!\n");
	return NULL;
  }

  header->magic           = be_int4(*(uint32_t *)(buffer+0));
  memcpy(header->version, buffer+4, 4);
  header->index_offset    = be_int8(*(uint64_t *)(buffer+8));
  header->index_len       = be_int4(*(uint32_t *)(buffer+16));
  header->nreads          = be_int4(*(uint32_t *)(buffer+20));
  header->header_len      = be_int2(*(uint16_t *)(buffer+24));
  header->key_len         = be_int2(*(uint16_t *)(buffer+26));
  header->flow_len        = be_int2(*(uint16_t *)(buffer+28));
  header->flowgram_format = be_int1(*(uint8_t *)(buffer+30));

  //Check file format and version
  if (header->magic != SFF_MAGIC || memcmp(header->version, SFF_VERSION, 4)) {
	free_header(header);
	return NULL;
  }

  // flow
  if(NULL == (header->flow = (char *)malloc((header->flow_len)+1))) {
    REprintf("Memory error!\n");
    free_header(header);
  } else {
    memcpy(header->flow, buffer+31, header->flow_len);
    header->flow[header->flow_len] = '\0';
  }

  // key
  if(NULL == (header->key = (char *)malloc((header->key_len)+1))) {
    REprintf("Memory error!\n");
    free_header(header);
  } else {
    memcpy(header->key, buffer+31+header->flow_len, header->key_len);
    header->key[header->key_len] = '\0';
  }

  //set pointer behind end of header
  pointer = header->header_len;


  /*****************************
   * Read in all read sections *
   *****************************/

  // Allocate space for container
  sff_read **reads;
  if(NULL == (reads = (sff_read **)calloc(header->nreads, sizeof(sff_read *)))) {
    REprintf("Memory error!\n");
    free_header(header);
    return NULL;
  }
  sff_container *container;
  if(NULL == (container = (sff_container *)calloc(1, sizeof(*container)))) {
    REprintf("Memory error!\n");
    free_header(header);
    free(reads);
	return NULL;
  }
  container->header = header;
  container->reads = reads;

  // Loop over all read sections until end of file
  int read_count = 0;
  while(1){

    if(pointer == header->index_offset) {
      //skip index
      pointer = pointer + header->index_len;
      continue;
    } else if(pointer >= fileLen){
      //end of file
      break;
    } else {
      /***********************
       * Read in read header *
       ***********************/
      sff_read *read = NULL;
      if(NULL == (read = (sff_read *)calloc(1, sizeof(*read)))) {
	    REprintf("Memory error!\n");
	    free_container(container, read_count);
	    return NULL;
	  }

      read->header_len         = be_int2(*(uint16_t *)(buffer+pointer+0));
      read->name_len           = be_int2(*(uint16_t *)(buffer+pointer+2));
      read->nbases             = be_int4(*(uint32_t *)(buffer+pointer+4));
      read->clip_qual_left     = be_int2(*(uint16_t *)(buffer+pointer+8));
      read->clip_qual_right    = be_int2(*(uint16_t *)(buffer+pointer+10));
      read->clip_adapter_left  = be_int2(*(uint16_t *)(buffer+pointer+12));
      read->clip_adapter_right = be_int2(*(uint16_t *)(buffer+pointer+14));

      // name
      if(NULL == (read->name  = (char *)malloc((read->name_len)+1))) {
        REprintf("Memory error!\n");
        free_container(container, read_count+1);
	    return NULL;
      } else {
	    memcpy(read->name, buffer+pointer+16, read->name_len);
	    read->name[read->name_len] = '\0';
	  }

	  //set pointer behind read header
	  pointer = pointer + read->header_len;

	  /*********************
	   * Read in read body *
	   *********************/
	  // flowgram
	  if(NULL == (read->flowgram = (uint16_t *)malloc(header->flow_len*2))) {
	    REprintf("Memory error!\n");
	    free_container(container, read_count+1);
	    return NULL;
      } else {
        memcpy(read->flowgram, buffer+pointer, header->flow_len*2);
        int k;
	    for(k=0;k<header->flow_len;k++) {
	      read->flowgram[k] = be_int2(read->flowgram[k]);
        }
	  }

      // flow index
      if(NULL == (read->flow_index = (uint8_t *)malloc(read->nbases))) {
        REprintf("Memory error!\n");
        free_container(container, read_count+1);
	    return NULL;
      } else {
        memcpy(read->flow_index, buffer+pointer+2*header->flow_len, read->nbases);
	  }

      // bases
      if(NULL == (read->bases = (char *)malloc((read->nbases)+1))) {
	    REprintf("Memory error!\n");
	    free_container(container, read_count+1);
	    return NULL;
      } else {
        memcpy(read->bases, buffer+pointer+2*header->flow_len+read->nbases, read->nbases);
        read->bases[read->nbases] = '\0';
	  }

      // quality
      if(NULL == (read->quality = (uint8_t *)malloc(read->nbases))) {
        REprintf("Memory error!\n");
        free_container(container, read_count+1);
	    return NULL;
      } else {
        memcpy(read->quality, buffer+pointer+2*header->flow_len+2*read->nbases, read->nbases);
	  }

	  //set pointer behind read and 8 byte padding
      pointer = pointer + 2*header->flow_len + 3*read->nbases;
      if(!(pointer % 8 == 0)) {
        pointer = pointer + (8 - (pointer % 8));
      }

      //add read into container
      container->reads[read_count] = read;

      read_count++;
    }
  }
  free(buffer);
  //printf("File successfully read! \n");
  return container;
}


/******************************************
 * Helper functions for memory management *
 ******************************************/

void free_header(sff_header *header) {
    if (!header)
	  return;
    if (header->flow)
	  free(header->flow);
    if (header->key)
	  free(header->key);
    free(header);
}

void free_read(sff_read *read) {
    if (!read)
	  return;
    if (read->name)
	  free(read->name);
    if (read->flowgram)
	  free(read->flowgram);
    if (read->flow_index)
	  free(read->flow_index);
    if (read->bases)
	  free(read->bases);
    if (read->quality)
	  free(read->quality);
    free(read);
}

void free_container(sff_container *container, int alloced_reads) {
  if (!container)
    return;
  if (container->header)
    free_header(container->header);
  if (container->reads) {
    int i;
    for(i=0; i<alloced_reads; i++) {
      if(container->reads[i]) {
        free_read(container->reads[i]);
      }
    }
    free(container->reads);
  }
  free(container);
}


/***************************************************************
 * Function returning a list in R-style given a .sff file name *
 ***************************************************************/

SEXP readSFFfromR(SEXP filename) {
  int nfiles;
  int i, j;
  SEXP flowgram_format, flow, key;
  SEXP r_names;
  SEXP clip_qual_left, clip_qual_right, clip_adapter_left, clip_adapter_right;
  SEXP flowgram_list, flow_index_list, quality_list, bases;
  SEXP list, list_names;

  PROTECT(filename = AS_CHARACTER(filename));
  nfiles = LENGTH(filename);
  char *filenames[nfiles];
  for(i = 0; i<nfiles; i++) {
    filenames[i] = R_alloc(strlen(CHAR(STRING_ELT(filename, i))), sizeof(char));
    strcpy(filenames[i], CHAR(STRING_ELT(filename, i)));
    //printf("file %d: %s\n", i, filenames[i]);
  }
  UNPROTECT(1);

  sff_container *co;
  co = readSFF(filenames[0]);
 //if NULL = co

  // flowgram_format
  PROTECT(flowgram_format = NEW_INTEGER(1));
  INTEGER(flowgram_format)[0] = co->header->flowgram_format;

  // flow
  PROTECT(flow = NEW_CHARACTER(1));
  SET_STRING_ELT(flow, 0, mkChar(co->header->flow));

  // key
  PROTECT(key = NEW_CHARACTER(1));
  SET_STRING_ELT(key, 0, mkChar(co->header->key));

  // character vector of read names
  PROTECT(r_names = NEW_CHARACTER(co->header->nreads));
  for(i=0; i<co->header->nreads; i++) {
    SET_STRING_ELT(r_names, i, mkChar(co->reads[i]->name));
  }

  // clipping
  PROTECT(clip_qual_left = NEW_INTEGER(co->header->nreads));
  for(i=0; i<co->header->nreads; i++) {
    INTEGER(clip_qual_left)[i] = co->reads[i]->clip_qual_left;
  }
  setAttrib(clip_qual_left, R_NamesSymbol, r_names);

  PROTECT(clip_qual_right = NEW_INTEGER(co->header->nreads));
  for(i=0; i<co->header->nreads; i++) {
    INTEGER(clip_qual_right)[i] = co->reads[i]->clip_qual_right;
  }
  setAttrib(clip_qual_right, R_NamesSymbol, r_names);

  PROTECT(clip_adapter_left = NEW_INTEGER(co->header->nreads));
  for(i=0; i<co->header->nreads; i++) {
    INTEGER(clip_adapter_left)[i] = co->reads[i]->clip_adapter_left;
  }
  setAttrib(clip_adapter_left, R_NamesSymbol, r_names);

  PROTECT(clip_adapter_right = NEW_INTEGER(co->header->nreads));
  for(i=0; i<co->header->nreads; i++) {
    INTEGER(clip_adapter_right)[i] = co->reads[i]->clip_adapter_right;
  }
  setAttrib(clip_adapter_right, R_NamesSymbol, r_names);

  // character vector of bases with read names as names
  PROTECT(bases = NEW_CHARACTER(co->header->nreads));
  for(i=0; i<co->header->nreads; i++) {
    SET_STRING_ELT(bases, i, mkChar(co->reads[i]->bases));
  }
  setAttrib(bases, R_NamesSymbol, r_names);

  // lists of flowgram, flow_index and quality with read names as names
  PROTECT(flowgram_list = allocVector(VECSXP, co->header->nreads));
  for(i = 0; i < co->header->nreads; i++) {
    SEXP flowgram;
    PROTECT(flowgram = NEW_INTEGER(co->header->flow_len));
    for(j=0; j<co->header->flow_len; j++) {
      INTEGER(flowgram)[j] = co->reads[i]->flowgram[j];
    }
    SET_VECTOR_ELT(flowgram_list, i, flowgram);
    UNPROTECT(1);
  }
  setAttrib(flowgram_list, R_NamesSymbol, r_names);

  PROTECT(flow_index_list = allocVector(VECSXP, co->header->nreads));
  for(i = 0; i < co->header->nreads; i++) {
    SEXP flow_index;
    PROTECT(flow_index = NEW_INTEGER(co->reads[i]->nbases));
    for(j=0; j<co->reads[i]->nbases; j++) {
      INTEGER(flow_index)[j] = co->reads[i]->flow_index[j];
    }
    SET_VECTOR_ELT(flow_index_list, i, flow_index);
    UNPROTECT(1);
  }
  setAttrib(flow_index_list, R_NamesSymbol, r_names);

  PROTECT(quality_list = allocVector(VECSXP, co->header->nreads));
  for(i = 0; i < co->header->nreads; i++) {
    SEXP quality;
    PROTECT(quality = NEW_INTEGER(co->reads[i]->nbases));
    for(j=0; j<co->reads[i]->nbases; j++) {
      INTEGER(quality)[j] = co->reads[i]->quality[j];
    }
    SET_VECTOR_ELT(quality_list, i, quality);
    UNPROTECT(1);
  }
  setAttrib(quality_list, R_NamesSymbol, r_names);

  // list representing the complete .sff file
  PROTECT(list = allocVector(VECSXP, 11));
  SET_VECTOR_ELT(list, 0, flowgram_format);
  SET_VECTOR_ELT(list, 1, flow);
  SET_VECTOR_ELT(list, 2, key);
  SET_VECTOR_ELT(list, 3, clip_qual_left);
  SET_VECTOR_ELT(list, 4, clip_qual_right);
  SET_VECTOR_ELT(list, 5, clip_adapter_left);
  SET_VECTOR_ELT(list, 6, clip_adapter_right);
  SET_VECTOR_ELT(list, 7, flowgram_list);
  SET_VECTOR_ELT(list, 8, flow_index_list);
  SET_VECTOR_ELT(list, 9, bases);
  SET_VECTOR_ELT(list, 10, quality_list);
  char *lnames[11] = {"flowgramFormat", "flowChars", "keySequence", "clipQualityLeft", "clipQualityRight", "clipAdapterLeft", "clipAdapterRight", "flowgrams", "flowIndexes", "reads", "qualityScores"};
  PROTECT(list_names = NEW_CHARACTER(11));
  for(i=0; i<11; i++) {
    SET_STRING_ELT(list_names, i, mkChar(lnames[i]));
  }
  setAttrib(list, R_NamesSymbol, list_names);
  UNPROTECT(14);

  return list;
}
