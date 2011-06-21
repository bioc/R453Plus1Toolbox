// sff_header holds all the information contained in the header of an .sff file
typedef struct {
    uint32_t  magic;
    char      version[4];
    uint64_t  index_offset;
    uint32_t  index_len;
    uint32_t  nreads;
    uint16_t  header_len;
    uint16_t  key_len;
    uint16_t  flow_len;
    uint8_t   flowgram_format;
    char     *flow;
    char     *key;
} sff_header;

// define the only values available for these two variable at the moment
#define SFF_MAGIC   0x2e736666 /* ".sff" */
#define SFF_VERSION "\0\0\0\1"

// sff_read holds the information for a single read, consisting of the header and the data section
typedef struct {
    uint16_t  header_len;
    uint16_t  name_len;
    uint32_t  nbases;
    uint16_t  clip_qual_left;
    uint16_t  clip_qual_right;
    uint16_t  clip_adapter_left;
    uint16_t  clip_adapter_right;
    char      *name;
    uint16_t  *flowgram;   // x 100.0 
    uint8_t   *flow_index; // relative to last
    char      *bases;
    uint8_t   *quality;
} sff_read;

// sff_container holds a single sff_header and multiple sff_reads and so represents a complete .sff file
typedef struct {
    sff_header  *header;
    sff_read    **reads;
} sff_container;

// Given the file name of a .sff file this function returns a sff_container
sff_container *readSFF(char *filename);

// Helper functions for memory management
void free_header(sff_header *header);
void free_read(sff_read *read);
void free_container(sff_container *container, int alloced_reads);
