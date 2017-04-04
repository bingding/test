/*
 * 
 * 
 * We need a function to make a fake bam1_t structure for unit-testing.
 * 
 * This structure requires the following components:
 * - CIGAR string (array of int32_t)
 * - Sequence string (2 nucleotides per byte)
 * - Quality string (1 quality value per byte)
 * - pos (int?)
 * - mate pos (int?)
 * - flag (int32_t)
 * - chr
 * - mate chr?
 * - read name
 * 
 * TASK:
 * 1) parse CIGAR string to BAM structure : 23M2S4D => array of 3 ints
 *      a) cut the string after each of the separators (MIDNSHP => 0123456)
 *      b) convert first part to int, and encode letter
 * 
 * MEMORY LAYOUT:
 * 
 * The cigar ops, the qname, and the sequence are stored in a single `data` array.
 * The layout is:
 * uint8_t * data;
 * sizeof(data) =  core->n_cigar * 4 (because 32 bits per operation) 
 *               + length of the qname from core->lqname (with \0)
 *               + length of the query sequence + 1  / 2 (2 nucleotides per byte)
 *               + length of the quality string (with \0)
 *               + aux (unused ?)
 * 
 * bam1_t
 * Structure for one alignment.
 *   typedef struct { 
 *       bam1_core_t core;  /// pointer to core, see below for members
 *       int l_aux, data_len, m_data; /// data_len is important!
 *       uint8_t *data;  // pointer to data
 *   } bam1_t;  
 *   
 *   Fields:
 *   core : core information about the alignment
 *   l_aux : length of auxiliary data
 *   data_len : current length of bam1_t::data
 *   m_data : maximum length of bam1_t::data
 *   data: all variable-length data, concatenated; structure: cigar-qname-seq-qual-aux 
 * 
 * 
 * bam1_core_t
 *
 *   Structure for core alignment information.
 *   typedef struct { 
 *       int32_t tid; // chromosome ID: need header!
 *       int32_t pos;  
 *       uint32_t bin:16, qual:8, l_qname:8; 
 *       uint32_t flag:16, n_cigar:16; 
 *       int32_t l_qseq; 
 *       int32_t mtid; //mate info
 *       int32_t mpos; // mate info
 *       int32_t isize; 
 *   } bam1_core_t;  
 *
 *   Fields:
 *  tid :chromosome ID, defined by bam_header_t
 *   pos:     0-based leftmost coordinate
 *   strand: strand; 0 for forward and 1 otherwise
 *   bin :     bin calculated by bam_reg2bin()
 *   qual :    mapping quality
 *   l_qname:     length of the query name
 *   flag:    bitwise flag
 *   n_cigar :    number of CIGAR operations
 *   l_qseq : length of the query sequence (read)
 * 
 * g++ -Wall helpers.cc -I /mnt/pipeline-programs/samtools/samtools-1.2/ -I /mnt/pipeline-programs/samtools/samtools-1.2/htslib-1.2.1/ /mnt/pipeline-programs/samtools/samtools-1.2/htslib-1.2.1/libhts.a -g -lz -lpthread
 */
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstdint>
#include "sam.h"
#include "htslib/sam.h"
#include "htslib/kstring.h"

namespace str_utils {
// Some split functions from stackoverflow
// https://stackoverflow.com/questions/236129/split-a-string-in-c
void split(const std::string &s, char delim, 
           std::vector<std::string> &elems,
           bool ignore_blank) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        if (ignore_blank && item.size() == 0) continue;
        elems.push_back(item);
    }
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems, true);
    return elems;
}

std::vector<std::string> split(const std::string & str, std::string delims)
{
    std::stringstream str_stream(str);
    std::string line;
    std::vector<std::string> out_arr;
    while(std::getline(str_stream, line)) 
    {
        std::size_t prev = 0, pos;
        while ((pos = line.find_first_of(delims, prev)) != std::string::npos)
        {
            if (pos > prev)
                out_arr.push_back(line.substr(prev, pos-prev));
            prev = pos+1;
        }
        if (prev < line.length())
            out_arr.push_back(line.substr(prev, std::string::npos));
    }
    return out_arr;
}


int chrom2number(std::string chr)
{
  int i;
  if(chr=="chrX")
    i=23;
  else if(chr=="chrY")
    i=24;
  else
    i=atoi((chr.substr(chr.find("chr")+3)).c_str());
  return(i-1);
}

}

namespace cigar_utils {
char num_to_cigar_op(int num) {
    return "MIDNSHP"[num];
}

int cigar_op_to_num(char op) 
{
    switch (op) {
        case 'M':
            return BAM_CMATCH;
            break;
        case 'I':
            return BAM_CINS;
            break;
        case 'D':
            return BAM_CDEL;
            break;
        case 'N':
            return BAM_CREF_SKIP;
            break;
        case 'S':
            return BAM_CSOFT_CLIP;
            break;
        case 'H':
            return BAM_CHARD_CLIP;
            break;
        case 'P':
            return BAM_CPAD;
            break;
        case '=':
        case 'X':
            std::cerr << "Not implemented CIGAR op: " << std::to_string(op) 
                      << std::endl;
            exit(1);
            break;
        default:
            std::cerr << "Unrecognized CIGAR op: " << std::to_string(op) 
                      << std::endl;
            exit(1);
            break;
    }
    return -1;
}



// Split string, from stackoverflow: https://stackoverflow.com/questions/7621727/split-a-string-into-words-by-multiple-delimiters
uint32_t split_cigar_string(std::string str, 
                            std::vector<std::pair<int, char>> & cigar_items) 
{
    std::stringstream string_stream(str);
    std::string line;
    const std::string delim = "MIDNSHP";
    uint32_t seq_length = 0;
    while(std::getline(string_stream, line)) 
    {
        std::size_t prev = 0, pos;
        while ((pos = line.find_first_of(delim, prev)) != std::string::npos)
        {
            if (pos > prev) {
                size_t converted_digits = 0;
                uint32_t op_count = std::stoi(line.substr(prev, pos-prev), 
                                         &converted_digits);
                char op_name = line.at(pos);
                if (converted_digits != pos - prev) {
                    std::cerr << "GIGAR string format error: unrecognized "
                                 "character at position " 
                              << pos - converted_digits << std::endl;
                    exit(1);
                }
                if (op_count > (UINT32_MAX >> 4)) {
                    std::cerr << "GIGAR string format error: number of items"
                                 "too big at pos " 
                              << pos << std::endl;
                }
                cigar_items.push_back(std::make_pair(op_count, op_name));
                seq_length += op_count;
            } else {
                std::cerr << "GIGAR string format error: delimiter without "
                             "length at pos " 
                          << pos << std::endl;
                exit(1);
            }
            prev = pos+1;
        }
        if (prev < line.length()) {
            std::string remaining = line.substr(prev, std::string::npos);
            std::cerr << "CIGAR string format error: extra characters after "
                         "last operation after pos " 
                      << prev << std::endl;
            exit(1);
        }
    }
    return seq_length;
}
}

namespace bam_utils {
/**
* @brief Compute a text-based sequence string from BAM data
* @param seq The sequence pointer
* @param length The sequence length
* 
* @returns A string containing the sequence
*/
std::string get_seq_str(uint8_t *seq, int length)
{
  std::stringstream line;
  //Get the sequence
  line.str(""); line.clear();
  for (int i=0;i<length;i++)
  {
    line<<bam_nt16_rev_table[bam1_seqi(seq,i)];
  }

  return line.str();
}

std::string get_qual_str(uint8_t * qseq, int length) 
{
    std::stringstream line;
    for (int i=0; i<length; i++) {
        line << qseq[i];
    }
    return line.str();
}
}



// Convert cigar items to an array of cigar values
int cigar_ops_to_array(const std::vector<std::pair<int, char>> & cigar_ops, 
                       uint32_t ** cigar_arr) 
{
    int n_cigar_ops = cigar_ops.size();
    (*cigar_arr) = (uint32_t *) malloc(n_cigar_ops * sizeof(uint32_t));
    for (int i = 0; i < n_cigar_ops; i++) {
        // 4 LSB are the CIGAR ops, other are length
        (*cigar_arr)[i] = (cigar_ops.at(i).first << BAM_CIGAR_SHIFT)
                          + cigar_utils::cigar_op_to_num(cigar_ops.at(i).second);
    }
    return n_cigar_ops;
}

// Just print to check
void print_cigar_array(uint32_t * cigar_arr, const int n_cigar) 
{
    for (int i=0; i < n_cigar; i++) {
        std::cout << (cigar_arr[i] >> BAM_CIGAR_SHIFT) << " " 
                  << cigar_utils::num_to_cigar_op(cigar_arr[i] & BAM_CIGAR_MASK) 
                  << std::endl;
    }
}


int seq_to_bam_seq(std::string seq, uint8_t ** bam_seq) 
{
    int n_items = (seq.size() + 1 ) / 2;
    int num = seq.size();
    (*bam_seq) = (uint8_t *) malloc(n_items * sizeof(uint8_t));
    for (int i=0; i<n_items; i++) {
        
        int nucl1;
        if (2*i < num) {
            nucl1 = bam_nt16_table[(uint8_t) seq.at(2 * i)];
        } else
            nucl1 = 0;
        int nucl2;
        if (2*i + 1 < num) {
            nucl2 = bam_nt16_table[(uint8_t) seq.at(2*i + 1)];
        }
        else
            nucl2 = 0;
        (*bam_seq)[i] = (nucl1 << 4) + nucl2;
    }
    return n_items;
}

// Warning: need exact sequence size
std::string bam_seq_to_seq(uint8_t * seq, const int seq_length) {
    std::string tmp;
    const int num_bytes = (seq_length + 1) / 2;
    for (int i=0; i < num_bytes; i++) {
        char curr = seq[i];
        uint8_t high = 0x0f & (curr >> BAM_CIGAR_SHIFT);
        uint8_t low = (curr & 0xf) | 0x00;
        tmp.push_back(bam_nt16_rev_table[high]);
        tmp.push_back(bam_nt16_rev_table[low]);
    }
    if (seq_length % 2 != 0) tmp.erase(tmp.length()-1);
    return tmp;
}


// Now make auxiliary fields
const std::string make_aux_string(std::string aux_name, std::string value)
{
    // Use The H or Z type for a normal string (ie barcode)
    const std::string aux_str = aux_name + 'Z' + value;
    return aux_str;
    
}

bam1_t * make_bam_record(std::string seq, std::string qual_seq,
                         std::string cigar_str, std::string q_name,
                         std::pair<std::string, std::string> aux_field,
                         int tid, int pos, int mtid, int mpos, int qual,
                         int flag
                        ) 
{
    if (seq.size() != qual_seq.size()) {
        std::cerr << "Invalid sequence or quality sequence: size of seq: " 
                  << seq.size() << " size of quality seq: " 
                  << qual_seq.size() << std::endl;
    }
    bam1_t * b = bam_init1();
    // make sequence string:
    uint8_t * b_seq = NULL;
    int b_seq_size = seq_to_bam_seq(seq, &b_seq);
    
    // make quality string
    const char * b_qseq = qual_seq.c_str();
    
    // make cigar array
    std::vector<std::pair<int, char>> cigar_items;
    uint32_t * cigar_array = NULL;
    auto seq_length_from_cigar = cigar_utils::split_cigar_string(cigar_str, 
                                                                 cigar_items);
    auto num_cigar_ops = cigar_ops_to_array(cigar_items, &cigar_array);
    
    if (seq_length_from_cigar != seq.size()) {
        std::cerr << "Invalid CIGAR or sequence: size of CIGAR: " 
                  << seq_length_from_cigar << " size of seq: " 
                  << seq.size() << std::endl;
    }
    
    // make qname
    const char * b_qname = q_name.c_str();
    
    // make aux field
    const std::string aux_str = make_aux_string(aux_field.first, aux_field.second);
    
    // Copy into data:
    // cigar_ops / qname / seq / qseq / aux
    const int data_length =   q_name.length() + 1 + num_cigar_ops * 4
                            + b_seq_size + qual_seq.size() + aux_str.size() + 1;
    
    b->data = (uint8_t *) malloc(data_length * sizeof(uint8_t));
    uint8_t * curr = b->data;
    memcpy(curr, b_qname, q_name.length() + 1);
    curr += q_name.length() + 1;
    memcpy(curr, cigar_array, num_cigar_ops * 4);
    curr += num_cigar_ops * 4;
    memcpy(curr, b_seq, b_seq_size);
    curr += b_seq_size;
    memcpy(curr, b_qseq, qual_seq.size());
    curr += qual_seq.size() ;
    memcpy(curr, aux_str.c_str(), aux_str.size() + 1 );
    
    // data length
    b->m_data = b->l_data = data_length;
    
    // Now populate other fields
    b->core.tid = tid;
    b->core.mtid = mtid;
    b->core.pos = pos;
    b->core.mpos = mpos;
    b->core.qual = qual;
    b->core.l_qname = q_name.size() + 1;
    b->core.n_cigar = num_cigar_ops;
    b->core.l_qseq = seq.size();
    
    // Unused:
    b->core.bin = 0;
    b->core.flag = 0; // FIXME TBD   
    b->core.isize = -1; // what is this?
    
    // cleanup seq string
    free(b_seq);
    // cleanup cigar array
    free(cigar_array);
    return b;
}


bam1_t * sam_read_to_bam(std::string sam_info) {
    // Cut sam read and push it to a strucuture
    //qname NB501251:16:HKGLWBGXX:1:11101:17614:1077 
    //flag 73
    //tid chr10
    //pos 123240469
    //mapq 37
    //cigar 72M
    //mate chr =
    //mate pos 123240469
    //insertsize 0
    //seq AAAATGAATAAAATTTATATACAATTATNTAAAATTCCAATTTCAGGATCCATAAAATTTTACTAGAAGACA
    //qseq EEEEEEEEEEEEEEEEEEEEEEEEEEEE#EEEEEEEEEAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEE/       
    //aux BC:Z:TTCCTnATATGACAGT   XT:A:U  NM:i:1  SM:i:37 AM:i:0  X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:28A43
    auto items = str_utils::split(sam_info, " \t");
    std::string qname = items.at(0);
    int flag = std::stoi(items.at(1));
    int tid = str_utils::chrom2number(items.at(2));
    int pos = std::stoi(items.at(3));
    int mapq = std::stoi(items.at(4));
    std::string cigar = items.at(5);
    int mtid = (items.at(6) == "=") ? tid : std::stoi(items.at(6));
    int mpos = std::stoi(items.at(7));
    int isize = 0; // FIXME unused
    std::string seq = items.at(9);
    std::string qseq = items.at(10);
    std::string bc;
    for (size_t i=11; i<items.size(); i++) {
        if (items.at(i).find("BC:Z") == 0) { // only check the barcode
            bc = items.at(i).substr(4);
        }
    }
    return make_bam_record(seq, qseq, cigar, qname, std::make_pair("BC", bc),
                         tid, pos, mtid, mpos, mapq, flag); 
    
}

std::vector<bam1_t*> sam_reads_to_bam(std::string reads)
{
    std::vector<bam1_t*> bam_reads;
    std::vector<std::string> sam_reads = str_utils::split(reads, "\n");
    for (const auto & read : sam_reads) {
        bam_reads.push_back(sam_read_to_bam(read));
    }
    return bam_reads;
}


int main(int argc, char * argv[]) 
{
    
    // Make/read CIGAR string
    std::vector<std::pair<int, char>> cigar_items;
    auto seq_length = cigar_utils::split_cigar_string("23M3D3I2S", cigar_items);
    for (const auto & item : cigar_items) {
        std::cout << item.first << " " << item.second << std::endl;
    }
    std::cout << "Sequence length should be " << seq_length << std::endl;
    uint32_t * cigar_array = NULL;
    auto num_cigar_ops = cigar_ops_to_array(cigar_items, &cigar_array);
    print_cigar_array(cigar_array, num_cigar_ops);
    free(cigar_array);
    
    //Make/read sequence string
    std::string seq = "ACTGA";
    uint8_t * b_seq = NULL;
    auto n_items = seq_to_bam_seq(seq, &b_seq);
    std::cout << "seq size : " << n_items << std::endl;
    std::cout << bam_seq_to_seq(b_seq, seq.size()) << std::endl;
    free(b_seq);
    
    // Make auxiliary field
    std::string aux_field = make_aux_string("BC", "TCACAGAG");
    
    /*
    bam1_t * my_record = make_bam_record(seq, "AAAAA", "2M3S", "myread",
                         std::make_pair("BC", "TCACAGAG"),
                         1, 1337, 2, 1400, 0xFF, 0x00);
    */
    /*
    bam1_t * my_record = sam_read_to_bam("NB501251:16:HKGLWBGXX:1:11101:19024:1078\t"
                                         "163\tchr7\t140484749\t29\t65M7S\t=\t140484816\t"
                                         "139\tACTTAAGTNTTTACAAATNAGTTNANCTANTNAATANAACATNGNAAGNTNANNTGGGAGANTGANATNNTA\t"
                                         "EEEEEEEE#EEEEEEEEE#EEEE#E#EEE#E#EEEE#EEEEE#E#EEE#E#E##EEAEEEE#AEE#EE##EE\t"
                                         "BC:Z:TCTTCnGTTGTATCCT\tXT:A:M\tNM:i:14\tSM:i:29\tAM:i:29\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:8T9C4T1A3C1A4A5G1C3A1G1A0A7A3");
    */
    bam1_t * my_record = sam_read_to_bam("NB501251:16:HKGLWBGXX:1:11101:17614:1077        73      chr10   123240469       37      72M     =       123240469       0       AAAATGAATAAAATTTATATACAATTATNTAAAATTCCAATTTCAGGATCCATAAAATTTTACTAGAAGACA  EEEEEEEEEEEEEEEEEEEEEEEEEEEE#EEEEEEEEEAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEE/        BC:Z:TTCCTnATATGACAGT   XT:A:U  NM:i:1  SM:i:37 AM:i:0  X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:28A43");
    
    
    // Now retrieve fields:
    std::cerr << "Barcode: " << bam_aux_get(my_record, "BC") << std::endl;
    std::cerr << "Read name: " << bam_get_qname(my_record) << std::endl;
    std::cerr << "Sequence: " << bam_utils::get_seq_str(bam_get_seq(my_record), 
                                             my_record->core.l_qseq) 
              << std::endl;
    std::cerr << "Quality Sequence: " << bam_utils::get_qual_str(bam_get_qual(my_record), 
                                             my_record->core.l_qseq)
              << std::endl;
    std::cerr << "data: " << my_record->data << std::endl;
    
    std::string reads = "NB501251:16:HKGLWBGXX:1:11101:17614:1077        73      chr10   123240469       37      72M     =       123240469       0       AAAATGAATAAAATTTATATACAATTATNTAAAATTCCAATTTCAGGATCCATAAAATTTTACTAGAAGACA  EEEEEEEEEEEEEEEEEEEEEEEEEEEE#EEEEEEEEEAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEE/        BC:Z:TTCCTnATATGACAGT   XT:A:U  NM:i:1  SM:i:37 AM:i:0  X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:28A43"
    "NB501251:16:HKGLWBGXX:1:11101:17614:1077        133     chr10   123240469       0       *       =       123240469       0       CGGCCCTGNNATCNCTGTNGCAANTNCGCNANTGTANCCTCANANTAGNCNCNNACAATATNTAANTANNTG  EEEEEEEE##EEE#EEEE#EEE6#E#EEE#A#EEEE#EEEEE#E#EEE#E#E##EEEEEEE#EEE#AE##EA        BC:Z:TTCCTnATATGACAGT"
    "NB501251:16:HKGLWBGXX:1:11101:20846:1077        89      chr10   123263353       37      72M     =       123263353       0       AGACCCCTGCCAGCATGGGGGTGTCTGCCGTTGAAGAGAGGCGNGTTGTTATCCTCACCAGCGGGGTGTTGG  /EEEEEEAEEEEEEEEEAEEEEEEEEEEEEEEEEEEEEEEEEE#EEEEEEEEEEEEEEEEEEEEEEEEEEEE        BC:Z:CTTGAnATTCTGACTT   XT:A:U  NM:i:1  SM:i:37 AM:i:0  X0:i:1  X1:i:0  XM:i:1  XO:i:0  XG:i:0  MD:Z:43T28"
    "NB501251:16:HKGLWBGXX:1:11101:20846:1077        181     chr10   123263353       0       *       =       123263353       0       NANNGGNATANTGGATTANNANTNTCTNTNTATGANCCTTNGNNACNANGGTCNTGGTNTGANNATGTCAGC  #/##A/#/E/#E6EEEE/##A#E#AEE#E#EAEEE#6EEE#A##//#A#EE/A#EEEE#AEE##AEE/EAEE        BC:Z:CTTGAnATTCTGACTT"
    "NB501251:16:HKGLWBGXX:1:11101:12836:1078        89      chr7    140495549       37      33M1D39M        =       140495549       0       AATTGATATGCAAATAATAAAATACTATGTACTAAAAAAAATTGGTATCGGCTCAATGAGAATTATAGGTAA  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEA        BC:Z:CATAGnCTCTGATATA   XT:A:U  NM:i:2  SM:i:37 AM:i:0  X0:i:1  X1:i:0  XM:i:1  XO:i:1  XG:i:1  MD:Z:24T8^A39"
    "NB501251:16:HKGLWBGXX:1:11101:12836:1078        181     chr7    140495549       0       *       =       140495549       0       AANNAANATTNTGAACTGNNGNGNTTCNGNTTCTGNCATGNTNTATNCNGAGANTGACNTTCANTGTTAGGC  AE##<A#EEE#EEEEEEE##E#E#E/E#E#EAEEE#EEEE#E#EEA#E#EEEE#EEEE#EEEE#EEEEEEEE        BC:Z:CATAGnCTCTGATATA"
    "NB501251:16:HKGLWBGXX:1:11101:22770:1078        73      chr4    1808514 37      72M     =       1808514 0       GCGGGAAGCGGCGGGGCTCACTCCTGAGCGCCCTGCCCGCAGGTACATGATCATGCGGGAGTGCTGGCATGC        EEEEEEEE/EEEEEEEEEEEEEEEEEAEAEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEAEEEEEEEA  BC:Z:TATAGnTGATCTGTGT   XT:A:U  NM:i:0  SM:i:37 AM:i:0  X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:72"
    "NB501251:16:HKGLWBGXX:1:11101:22770:1078        133     chr4    1808514 0       *       =       1808514 0       AGCACTCANGTCGNTGGANGTCANGNTAANGNCACGNTCCAGNTNCTCNANCNNCTGCTTGNAGGNGGNNCT        EEEEEEEE#EEEE#EEEE#AEEE#A#EEE#E#EEEE#EEEEE#E#EEE#A#A##EEEEEEE#EEE#EE##EE  BC:Z:TATAGnTGATCTGTGT"
    "NB501251:16:HKGLWBGXX:1:11101:19024:1078        83      chr7    140484816       29      72M     =       140484749       -139    TTATACTGTAAACTTTCTTCTATTAGGACATAATTAAAAGGATAAATTCTCAAATTTCTGAAGTACACTGAA  EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE        BC:Z:TCTTCnGTTGTATCCT   XT:A:U  NM:i:0  SM:i:29 AM:i:29 X0:i:1  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:72"
    "NB501251:16:HKGLWBGXX:1:11101:19024:1078        163     chr7    140484749       29      65M7S   =       140484816       139     ACTTAAGTNTTTACAAATNAGTTNANCTANTNAATANAACATNGNAAGNTNANNTGGGAGANTGANATNNTA  EEEEEEEE#EEEEEEEEE#EEEE#E#EEE#E#EEEE#EEEEE#E#EEE#E#E##EEAEEEE#AEE#EE##EE        BC:Z:TCTTCnGTTGTATCCT   XT:A:M  NM:i:14 SM:i:29 AM:i:29 XM:i:0  XO:i:0  XG:i:0  MD:Z:8T9C4T1A3C1A4A5G1C3A1G1A0A7A3";
    std::vector<bam1_t*> bb = sam_reads_to_bam(reads);
    
    for (size_t i=0; i<bb.size(); i++) {
        bam_destroy1(bb.at(i));
    }
    bam_destroy1(my_record);
    
    
    std::string bam_hdr = 
"@HD	VN:1.3	SO:coordinate\n"
"@SQ	SN:chr1	LN:249250621\n"
"@SQ	SN:chr2	LN:243199373\n"
"@SQ	SN:chr3	LN:198022430\n"
"@SQ	SN:chr4	LN:191154276\n"
"@SQ	SN:chr5	LN:180915260\n"
"@SQ	SN:chr6	LN:171115067\n"
"@SQ	SN:chr7	LN:159138663\n"
"@SQ	SN:chr8	LN:146364022\n"
"@SQ	SN:chr9	LN:141213431\n"
"@SQ	SN:chr10	LN:135534747\n"
"@SQ	SN:chr11	LN:135006516\n"
"@SQ	SN:chr12	LN:133851895\n"
"@SQ	SN:chr13	LN:115169878\n"
"@SQ	SN:chr14	LN:107349540\n"
"@SQ	SN:chr15	LN:102531392\n"
"@SQ	SN:chr16	LN:90354753\n"
"@SQ	SN:chr17	LN:81195210\n"
"@SQ	SN:chr18	LN:78077248\n"
"@SQ	SN:chr19	LN:59128983\n"
"@SQ	SN:chr20	LN:63025520\n"
"@SQ	SN:chr21	LN:48129895\n"
"@SQ	SN:chr22	LN:51304566\n"
"@SQ	SN:chrX	LN:155270560\n"
"@SQ	SN:chrY	LN:59373566\n"
"@PG	ID:bwa	PN:bwa	VN:0.6.1-r104-tpx\n"
"@PG	ID:bwa-715317FF	PN:bwa	VN:0.6.1-r104-tpx\n";
    
    bam_hdr_t * hdr = sam_hdr_parse(bam_hdr.size(), bam_hdr.c_str());
    std::cerr << "Number of targets: " << hdr->n_targets<< std::endl;
    std::string srd = "NB500878:160:HF72JAFXX:1:11312:7296:11299\t161\tchr1\t10034\t0\t80M\tchr10\t123180036\t0\tCCCTAACCCTAACCCTAACCCTTACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTACCCCAGCCC\tAAAAA6EEEEA6EEEEAE/<EEEEEEEEE//EEE//EEE///EE6/</6<E/</EEE////6/////<E////<E///A/\tXT:A:R\tNM:i:3\tSM:i:0\tAM:i:0\tX0:i:3\tX1:i:6\tXM:i:3\tXO:i:0\tXG:i:0\tMD:Z:22A48A4A3\tXA:Z:chr12,-95584,80M,3;chr15,-102521285,80M,3;chr4,-191043973,80M,4;chr4,-191044156,80M,4;chr12,-95452,80M,4;chr2,+243152505,80M,4;chr12,-95458,80M,4;chr1,+10040,80M,4;\tBC:Z:TTCCTnATATGACAGT\n";
    kstring_t s = {0, 0 , NULL};
    kputs(srd.c_str(), &s);
    bam1_t * bs = bam_init1();
    int res = sam_parse1(&s, hdr, bs);
    std::cerr << res << std::endl;
    std::cerr << "NEW: barcode: " << bam_aux_get(bs, "BC") << std::endl;
    return 0;
}
