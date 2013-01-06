#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "io_lib/os.h"
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "writeSFF.h"

/*******************************************************************
 * Helper function: get the list element named str, or return NULL *
********************************************************************/
SEXP getListElement(SEXP list, const char *str)
{
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

  for (R_len_t i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

/*******************************************
 * Function to write a file in .sff format *
********************************************/

SEXP writeSFFfromR(SEXP ll, SEXP fname)
{
  char *filename;
  FILE *file;
  PROTECT(fname = AS_CHARACTER(fname));
  filename = R_alloc(strlen(CHAR(STRING_ELT(fname, 0))), sizeof(char));
  strcpy(filename, CHAR(STRING_ELT(fname, 0)));
  file = fopen(filename, "wb+");
  
  PROTECT(ll = AS_LIST(ll));
  
  /* write header section */
  // magic number
  uint32_t magicNumber;
  magicNumber = be_int4(SFF_MAGIC);
  fwrite(&magicNumber, 4, 1, file);
  
  // sff version
  fwrite(SFF_VERSION, 4, 1, file);
  
  // index offset and length, no yet used
  uint64_t indexOffset = 0;
  fwrite(&indexOffset, 8, 1, file);
  uint32_t indexLength = 0;
  fwrite(&indexLength, 4, 1, file);
  
  // number of reads
  SEXP reads;
  PROTECT(reads = AS_CHARACTER(getListElement(ll, "reads")));
  uint32_t nreads = be_int4(GET_LENGTH(reads));
  fwrite(&nreads, 4, 1, file);

  // key length + key sequence
  SEXP kSequence;
  PROTECT(kSequence = getListElement(ll, "keySequence"));
  int kLength = strlen(CHAR(STRING_ELT(kSequence, 0)));
  char* keySequence = R_alloc(kLength, sizeof(char));
  strcpy(keySequence, CHAR(STRING_ELT(kSequence, 0)));
  uint16_t keyLength = be_int2(kLength);
  
  // number of flows, flowgram format and flow chars
  SEXP fChars;
  PROTECT(fChars = getListElement(ll, "flowChars"));
  int nFlows = strlen(CHAR(STRING_ELT(fChars, 0)));
  char* flowChars = R_alloc(nFlows, sizeof(char));
  strcpy(flowChars, CHAR(STRING_ELT(fChars, 0)));
  uint16_t numberFlows = be_int2(nFlows); 
  
  // headerlength
  int hlength = 31 + kLength + nFlows;
  int hlengthmod = hlength % 8;
  if(hlengthmod != 0)
  {
    hlength = hlength + (8 - hlengthmod);
  }
  uint16_t headerlength = be_int2(hlength);
	fwrite(&headerlength, 2, 1, file);
  
  fwrite(&keyLength, 2, 1, file);
  fwrite(&numberFlows, 2, 1, file);
  
  uint8_t flowgramFormat = be_int1(INTEGER(getListElement(ll, "flowgramFormat"))[0]);
  fwrite(&flowgramFormat, 1, 1, file);
  
  fwrite(flowChars, nFlows, 1, file);
  fwrite(keySequence, kLength, 1, file);
  
  // eight byte padding
  uint64_t filling = 0;
  fwrite(&filling, hlength - 31 - kLength - nFlows, 1, file);
  
  /* write read sections */
  // read header variables
  SEXP cql, cqr, cal, car;
  PROTECT(cql = AS_INTEGER(getListElement(ll, "clipQualityLeft")));
  PROTECT(cqr = AS_INTEGER(getListElement(ll, "clipQualityRight")));
  PROTECT(cal = AS_INTEGER(getListElement(ll, "clipAdapterLeft")));
  PROTECT(car = AS_INTEGER(getListElement(ll, "clipAdapterRight")));
  SEXP names = getAttrib(reads, R_NamesSymbol);
  char* rname;
  
  // read data variables
  SEXP fgrams, fIndexes, qScores;
  PROTECT(fgrams = AS_LIST(getListElement(ll, "flowgrams")));
  PROTECT(fIndexes = AS_LIST(getListElement(ll, "flowIndexes")));
  PROTECT(qScores = AS_LIST(getListElement(ll, "qualityScores"))); 
  int* flowgram;
  int* flowIndexes;
  int* qualityScores;
  
  for(int i = 0; i < GET_LENGTH(reads); ++i) {
  
    // read header length
    int rhlength = 16 + strlen(CHAR(STRING_ELT(names, i)));
    int rhlengthmod = rhlength % 8;
    if(rhlengthmod != 0)
    {
      rhlength = rhlength + (8 - rhlengthmod);
    }
    uint16_t readHeaderLength = be_int2(rhlength);
    fwrite(&readHeaderLength, 2, 1, file);
    
    //name length
    uint16_t nameLength = be_int2(strlen(CHAR(STRING_ELT(names, i))));
    fwrite(&nameLength, 2, 1, file);
    
    // number bases
    uint32_t numberBases = be_int4(strlen(CHAR(STRING_ELT(reads, i))));
    fwrite(&numberBases, 4, 1, file);
    
    //clipping quality and adapter
    uint16_t clipQualityLeft = be_int2(INTEGER(cql)[i]);
    fwrite(&clipQualityLeft, 2, 1, file);
    uint16_t clipQualityRight = be_int2(INTEGER(cqr)[i]);
    fwrite(&clipQualityRight, 2, 1, file);
    uint16_t clipAdapterLeft = be_int2(INTEGER(cal)[i]);
    fwrite(&clipAdapterLeft, 2, 1, file);
    uint16_t clipAdapterRight = be_int2(INTEGER(car)[i]);
    fwrite(&clipAdapterRight, 2, 1, file);
    
    //read name
    rname = R_alloc(strlen(CHAR(STRING_ELT(names, i))), sizeof(char));
    strcpy(rname, CHAR(STRING_ELT(names, i)));
    fwrite(rname, strlen(CHAR(STRING_ELT(names, i))), 1, file);
    
    //padding
    uint64_t rh_filling = 0;
    fwrite(&rh_filling, rhlength - 16 - strlen(CHAR(STRING_ELT(names, i))), 1, file);
    
    
    //read data section
    
    //flowgram values
    flowgram = INTEGER(VECTOR_ELT(fgrams, i));
    for(int j = 0; j < nFlows; j++) {
      uint16_t fg = be_int2(flowgram[j]);
      fwrite(&fg, 2, 1, file);
  	}

    //flowIndexes
    flowIndexes = INTEGER(VECTOR_ELT(fIndexes, i));
    for(int k = 0; k < strlen(CHAR(STRING_ELT(reads, i))); k++) {
      uint8_t fi = be_int1(flowIndexes[k]);
      fwrite(&fi, 1, 1, file);
    }
  
    //bases
    char* bases;
    bases = R_alloc(strlen(CHAR(STRING_ELT(reads, i))), sizeof(char));
    strcpy(bases, CHAR(STRING_ELT(reads, i)));
    fwrite(bases, strlen(CHAR(STRING_ELT(reads, i))), 1, file);
  
    //qualityScores
    qualityScores = INTEGER(VECTOR_ELT(qScores, i));
    for(int l = 0; l < strlen(CHAR(STRING_ELT(reads, i))); l++) {
      uint8_t qs = be_int1(qualityScores[l]);
      fwrite(&qs, 1, 1, file);
    }
  
    //padding
    uint64_t rd_filling = 0;
  	fwrite(&rd_filling, 8 - ((3 * strlen(CHAR(STRING_ELT(reads, i))) + 2 * nFlows) % 8), 1, file);
  }
  
  UNPROTECT(12);
  
  fclose(file);
  return(R_NilValue);
}
