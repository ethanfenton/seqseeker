/*
 * kmer_search.c — Fast k-mer–based sequence detector for paired gzipped FASTQ files
 *
 * Compile:
 *   gcc -O3 -march=native -o kmer_search kmer_search.c -lz
 *
 * Usage:
 *   kmer_search -1 R1.fastq.gz -2 R2.fastq.gz -o out.csv \
 *               -q SEQUENCE|-f seq.fasta \
 *               [-s sample] [-l lane] [-k 20] [-b 16] [-u 12] [-m 1] [-r 1000000]
 *               [-R full|5prime|3prime] [-L region_len]
 *
 * Output CSV columns:
 *   sample, lane, barcode, umi, n_kmer_hits, target_pos, strand, read_name, r2_seq
 *
 * The search sequence can be provided three ways:
 *   -q ACGT...    raw nucleotide string on the command line
 *   -f FILE       path to a FASTA file (first sequence record is used)
 *   -n NAME       built-in sequence name (see below)
 *
 * Region options control which portion of the target is indexed:
 *   -R full       entire sequence (default)
 *   -R 5prime     first -L bp (default 200)
 *   -R 3prime     last  -L bp (default 200)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <time.h>
#include <sys/stat.h>
#include <zlib.h>

/* ── Built-in sequences (macros so they can appear in struct initializers) ── */

#define SEQ_EGFP \
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAAC" \
    "GGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCT" \
    "GCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGC" \
    "CGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCAC" \
    "CATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACC" \
    "GCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAAC" \
    "AGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACAT" \
    "CGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGC" \
    "CCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTG" \
    "CTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA"

/* Venus/YFP — key mutations T203Y,S65G,V68L,S72A on eGFP backbone */
#define SEQ_VENUS \
    "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAAC" \
    "GGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCT" \
    "GCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGC" \
    "CGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCAC" \
    "CATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACC" \
    "GCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAAC" \
    "AGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACAT" \
    "CGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGC" \
    "CCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTG" \
    "CTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAT"

#define SEQ_MCHERRY \
    "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGG" \
    "GCAGCATGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCC" \
    "AAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCC" \
    "AAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAG" \
    "CGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATC" \
    "TACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAG" \
    "GCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGAC" \
    "GGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAAC" \
    "GTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGC" \
    "CGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"

#define SEQ_TDTOMATO \
    "ATGGTGAGCAAGGGCGAGGAGGTCATCAAAGAGTTCATGCGCTTCAAGGTGCGCATGGAAGGCTCCATGAACGG" \
    "CCACGAGTTTGAGATCGAAGGCGAAGGCGAAGGCCGCCCCTACGAAGGCACCCAGACCGCCAAGCTGAAGGTGAC" \
    "CAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCGCCTCAGTTCATGTACGGCTCCAGGGCCTTCGTGAAG" \
    "CACCCCGCCGACATCCCCGACTACAAGAAGCTGTCCTTCCCCGAAGGCTTCAAGTGGGAGCGCGTGATGAACTTCG" \
    "AGGATGGCGGCGTGGTGACCGTGACACAGGACTCCTCCTTACAAGATGGCGAATTCATCTACAAGGTTAAGCTTCG" \
    "CGGCACCAACTTCCCATCCGATGGCCCTGTAATGCAAAAGAAAACAATGGGCTGGGAAGCCTCCACCGAGCGCATG" \
    "TACCCTGAGGATGGCGCCTTGAAGGGCGAAATAAGCCAGAGGTTGAAACTGAAGGACGGCGGTCACTACGATGCCG" \
    "AGGTCAAGACCACGTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAATATCAAGTTGGACATC" \
    "ACCTCCCATAACGAGGATTATACCATCGTTGAGCAGTACGAGCGGGCCGAAGGCCGGCACTCCACCGGCGGCATGG" \
    "ACGAGCTGTACAAGGCGGCCATGGTGAGCAAGGGCGAGGAGGTCATCAAAGAGTTCATGCGCTTCAAGGTGCGCAT" \
    "GGAAGGCTCCATGAACGGCCACGAGTTTGAGATCGAAGGCGAAGGCGAAGGCCGCCCCTACGAAGGCACCCAGACC" \
    "GCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCGCCTCAGTTCATGTACGGCT" \
    "CCAGGGCCTTCGTGAAGCACCCCGCCGACATCCCCGACTACAAGAAGCTGTCCTTCCCCGAAGGCTTCAAGTGGGAG" \
    "CGCGTGATGAACTTCGAGGATGGCGGCGTGGTGACCGTGACACAGGACTCCTCCTTACAAGATGGCGAATTCATCTA" \
    "CAAGGTTAAGCTTCGCGGCACCAACTTCCCATCCGATGGCCCTGTAATGCAAAAGAAAACAATGGGCTGGGAAGCCT" \
    "CCACCGAGCGCATGTACCCTGAGGATGGCGCCTTGAAGGGCGAAATAAGCCAGAGGTTGAAACTGAAGGACGGCGGT" \
    "CACTACGATGCCGAGGTCAAGACCACGTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAATAT" \
    "CAAGTTGGACATCACCTCCCATAACGAGGATTATACC" \
    "ATCGTTGAGCAGTACGAGCGGGCCGAAGGCCGGCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA"

typedef struct { const char *name; const char *seq; } BuiltinEntry;
static const BuiltinEntry BUILTINS[] = {
    {"egfp",     SEQ_EGFP},
    {"venus",    SEQ_VENUS},
    {"mcherry",  SEQ_MCHERRY},
    {"tdtomato", SEQ_TDTOMATO},
    {NULL, NULL}
};

static const char *lookup_builtin(const char *name)
{
    for (int i = 0; BUILTINS[i].name; i++)
        if (strcasecmp(BUILTINS[i].name, name) == 0)
            return BUILTINS[i].seq;
    return NULL;
}

/* ── Constants ───────────────────────────────────────────────────────────── */
#define MAX_SEQ_LEN   100000
#define BUF_SIZE      8192
#define HASH_BITS     17        /* 131072 slots — comfortable for 100 kb targets */
#define HASH_SIZE     (1 << HASH_BITS)
#define HASH_MASK     (HASH_SIZE - 1)
#define GZBUF_BYTES   (1 << 21)

typedef enum { REGION_FULL = 0, REGION_5PRIME, REGION_3PRIME } Region;

/* ── Hash-table entry ─────────────────────────────────────────────────────── */
typedef struct {
    uint64_t kmer;
    int32_t  target_pos;
    int8_t   strand;
    uint8_t  occupied;
} Slot;

static Slot *g_htable = NULL;   /* heap-allocated in main */

/* ── Base helpers ─────────────────────────────────────────────────────────── */
static inline int base_to_bit(char c)
{
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default:             return -1;
    }
}

static inline char comp(char c)
{
    switch (c) {
        case 'A': return 'T'; case 'a': return 't';
        case 'C': return 'G'; case 'c': return 'g';
        case 'G': return 'C'; case 'g': return 'c';
        case 'T': return 'A'; case 't': return 'a';
        default:  return 'N';
    }
}

static void rev_comp(const char *src, char *dst, int len)
{
    for (int i = 0; i < len; i++)
        dst[len - 1 - i] = comp(src[i]);
    dst[len] = '\0';
}

/* ── Hash function ──────────────────────────────────────────────────────── */
static inline uint32_t kmer_hash(uint64_t v)
{
    v ^= v >> 33;
    v *= UINT64_C(0xff51afd7ed558ccd);
    v ^= v >> 33;
    v *= UINT64_C(0xc4ceb9fe1a85ec53);
    v ^= v >> 33;
    return (uint32_t)(v & HASH_MASK);
}

static void htable_insert(uint64_t kmer, int32_t pos, int8_t strand)
{
    uint32_t h = kmer_hash(kmer);
    while (g_htable[h].occupied && g_htable[h].kmer != kmer)
        h = (h + 1) & HASH_MASK;
    if (!g_htable[h].occupied) {
        g_htable[h].kmer       = kmer;
        g_htable[h].target_pos = pos;
        g_htable[h].strand     = strand;
        g_htable[h].occupied   = 1;
    }
}

static const Slot *htable_lookup(uint64_t kmer)
{
    uint32_t h = kmer_hash(kmer);
    while (g_htable[h].occupied) {
        if (g_htable[h].kmer == kmer) return &g_htable[h];
        h = (h + 1) & HASH_MASK;
    }
    return NULL;
}

/* ── Build k-mer index ───────────────────────────────────────────────────── */
static int build_kmer_hash(const char *seq, int seq_len, int k,
                           Region region, int region_len)
{
    if (k < 1 || k > 31) {
        fprintf(stderr, "k must be between 1 and 31\n");
        return -1;
    }

    /* Determine the slice of the target sequence to index */
    int start = 0, end = seq_len;
    if (region == REGION_5PRIME) {
        end = (region_len > 0 && region_len < seq_len) ? region_len : seq_len;
    } else if (region == REGION_3PRIME) {
        start = (region_len > 0 && seq_len - region_len > 0)
                ? seq_len - region_len : 0;
    }
    int slice_len = end - start;

    if (slice_len < k) {
        fprintf(stderr, "[WARN] Region length (%d) < k (%d); "
                        "no k-mers will be indexed. Widening to full.\n",
                slice_len, k);
        start = 0; end = seq_len; slice_len = seq_len;
    }

    uint64_t mask = (k < 32) ? (((uint64_t)1 << (2*k)) - 1) : UINT64_MAX;

    char *rc = malloc(seq_len + 1);
    if (!rc) { fprintf(stderr, "OOM\n"); return -1; }
    rev_comp(seq, rc, seq_len);

    int n = 0;

    /* Forward strand over [start, end) */
    for (int i = start; i <= end - k; i++) {
        uint64_t val = 0; int ok = 1;
        for (int j = 0; j < k; j++) {
            int b = base_to_bit(seq[i + j]);
            if (b < 0) { ok = 0; break; }
            val = ((val << 2) | (uint64_t)b) & mask;
        }
        if (ok) { htable_insert(val, i, +1); n++; }
    }

    /* Reverse-complement strand: rc[i] corresponds to fwd pos seq_len-i-k */
    for (int i = (seq_len - end); i <= (seq_len - start) - k; i++) {
        uint64_t val = 0; int ok = 1;
        for (int j = 0; j < k; j++) {
            int b = base_to_bit(rc[i + j]);
            if (b < 0) { ok = 0; break; }
            val = ((val << 2) | (uint64_t)b) & mask;
        }
        if (ok) { htable_insert(val, seq_len - i - k, -1); n++; }
    }

    free(rc);
    return n;
}

/* ── Sequence loading ─────────────────────────────────────────────────────── */

/* Strip non-ACGTN characters (spaces, newlines, numbers from FASTA) */
static int clean_sequence(const char *raw, char *out, int max_len)
{
    int n = 0;
    for (int i = 0; raw[i] && n < max_len; i++) {
        char c = toupper((unsigned char)raw[i]);
        if (c=='A'||c=='C'||c=='G'||c=='T'||c=='N') out[n++] = c;
    }
    out[n] = '\0';
    return n;
}

/* Load first sequence from a FASTA file; returns length or -1 on error */
static int load_fasta(const char *path, char *out, int max_len)
{
    FILE *fh = fopen(path, "r");
    if (!fh) { fprintf(stderr, "Cannot open FASTA: %s\n", path); return -1; }

    char line[4096];
    int in_seq = 0, n = 0;
    while (fgets(line, sizeof(line), fh)) {
        if (line[0] == '>') {
            if (in_seq) break;   /* stop after first record */
            in_seq = 1;
            continue;
        }
        if (!in_seq) continue;
        for (int i = 0; line[i] && n < max_len; i++) {
            char c = toupper((unsigned char)line[i]);
            if (c=='A'||c=='C'||c=='G'||c=='T'||c=='N') out[n++] = c;
        }
    }
    out[n] = '\0';
    fclose(fh);
    return n;
}

/* ── Read scan ───────────────────────────────────────────────────────────── */
static int scan_read(const char *seq, int slen, int k, uint64_t mask,
                     int *first_pos, int *first_strand)
{
    int hits = 0;
    uint64_t val = 0;
    int fill = 0;
    *first_pos = -1; *first_strand = 0;

    for (int i = 0; i < slen; i++) {
        int b = base_to_bit(seq[i]);
        if (b < 0) { val = 0; fill = 0; continue; }
        val = ((val << 2) | (uint64_t)b) & mask;
        if (++fill >= k) {
            const Slot *s = htable_lookup(val);
            if (s) {
                hits++;
                if (hits == 1) { *first_pos = s->target_pos; *first_strand = s->strand; }
            }
        }
    }
    return hits;
}

/* ── FASTQ reading ──────────────────────────────────────────────────────── */
static int read_fastq_record(gzFile fh,
                              char *header, char *seq,
                              char *plus,   char *qual)
{
    if (!gzgets(fh, header, BUF_SIZE)) return 0;
    if (!gzgets(fh, seq,    BUF_SIZE)) return 0;
    if (!gzgets(fh, plus,   BUF_SIZE)) return 0;
    if (!gzgets(fh, qual,   BUF_SIZE)) return 0;
    seq[strcspn(seq, "\r\n")] = '\0';
    return 1;
}

/* ── Helpers ─────────────────────────────────────────────────────────────── */
static void fmt_time(double sec, char *buf, int n)
{
    if (sec < 0 || sec > 604800.0) { snprintf(buf, n, "??:??:??"); return; }
    snprintf(buf, n, "%02d:%02d:%02d",
             (int)(sec/3600), ((int)(sec/60))%60, ((int)sec)%60);
}

static long get_file_size(const char *path)
{
    struct stat st;
    return (stat(path, &st) == 0) ? (long)st.st_size : -1L;
}

static void csv_field(FILE *out, const char *s, int n)
{
    for (int i = 0; s[i] && i < n; i++) {
        char c = s[i];
        fputc((c==',' || c=='\n' || c=='\r') ? '_' : c, out);
    }
}

static void print_builtins(void)
{
    fprintf(stderr, "Built-in sequence names (-n):\n");
    for (int i = 0; BUILTINS[i].name; i++)
        fprintf(stderr, "  %s  (%d bp)\n",
                BUILTINS[i].name, (int)strlen(BUILTINS[i].seq));
}

/* ═══════════════════════════════════════════════════════════════════════════ */
int main(int argc, char **argv)
{
    const char *r1_path         = NULL;
    const char *r2_path         = NULL;
    const char *sample          = "unknown";
    const char *lane            = "unknown";
    const char *outpath         = NULL;
    const char *seq_raw         = NULL;   /* -q raw sequence string */
    const char *seq_file        = NULL;   /* -f FASTA file          */
    const char *seq_name        = NULL;   /* -n built-in name       */
    int         kmer_size       = 20;
    int         barcode_len     = 16;
    int         umi_len         = 12;
    int         min_kmers       = 1;
    long        report_interval = 1000000L;
    Region      region          = REGION_FULL;
    int         region_len      = 200;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i],"-1") && i+1<argc) r1_path         = argv[++i];
        else if (!strcmp(argv[i],"-2") && i+1<argc) r2_path         = argv[++i];
        else if (!strcmp(argv[i],"-s") && i+1<argc) sample          = argv[++i];
        else if (!strcmp(argv[i],"-l") && i+1<argc) lane            = argv[++i];
        else if (!strcmp(argv[i],"-o") && i+1<argc) outpath         = argv[++i];
        else if (!strcmp(argv[i],"-q") && i+1<argc) seq_raw         = argv[++i];
        else if (!strcmp(argv[i],"-f") && i+1<argc) seq_file        = argv[++i];
        else if (!strcmp(argv[i],"-n") && i+1<argc) seq_name        = argv[++i];
        else if (!strcmp(argv[i],"-k") && i+1<argc) kmer_size       = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-b") && i+1<argc) barcode_len     = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-u") && i+1<argc) umi_len         = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-m") && i+1<argc) min_kmers       = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-r") && i+1<argc) report_interval = atol(argv[++i]);
        else if (!strcmp(argv[i],"-L") && i+1<argc) region_len      = atoi(argv[++i]);
        else if (!strcmp(argv[i],"-R") && i+1<argc) {
            const char *rv = argv[++i];
            if      (!strcmp(rv,"5prime"))  region = REGION_5PRIME;
            else if (!strcmp(rv,"3prime"))  region = REGION_3PRIME;
            else if (!strcmp(rv,"full"))    region = REGION_FULL;
            else { fprintf(stderr,"Unknown region: %s\n",rv); return 1; }
        }
        else if (!strcmp(argv[i],"--list-sequences")) { print_builtins(); return 0; }
        else { fprintf(stderr, "Unknown argument: %s\n", argv[i]); return 1; }
    }

    if (!r1_path || !r2_path || !outpath) {
        fprintf(stderr,
            "Usage: %s -1 R1.fastq.gz -2 R2.fastq.gz -o out.csv\n"
            "       ( -q SEQUENCE | -f seq.fasta | -n NAME )\n"
            "       [-s sample] [-l lane] [-k %d] [-b %d] [-u %d]\n"
            "       [-m %d] [-r %ld] [-R full|5prime|3prime] [-L %d]\n"
            "       [--list-sequences]\n",
            argv[0], kmer_size, barcode_len, umi_len,
            min_kmers, report_interval, region_len);
        return 1;
    }

    /* Load target sequence */
    static char target_seq[MAX_SEQ_LEN + 1];
    int target_len = 0;

    if (seq_name) {
        const char *b = lookup_builtin(seq_name);
        if (!b) {
            fprintf(stderr, "Unknown built-in sequence: %s\n", seq_name);
            print_builtins();
            return 1;
        }
        target_len = clean_sequence(b, target_seq, MAX_SEQ_LEN);
    } else if (seq_file) {
        target_len = load_fasta(seq_file, target_seq, MAX_SEQ_LEN);
        if (target_len <= 0) return 1;
    } else if (seq_raw) {
        target_len = clean_sequence(seq_raw, target_seq, MAX_SEQ_LEN);
    } else {
        fprintf(stderr, "No target sequence provided. Use -q, -f, or -n.\n");
        return 1;
    }

    if (target_len < kmer_size) {
        fprintf(stderr, "Target sequence (%d bp) shorter than k (%d).\n",
                target_len, kmer_size);
        return 1;
    }

    fprintf(stderr, "[%s/%s] Target: %d bp  region: %s",
            sample, lane, target_len,
            region==REGION_FULL ? "full" :
            region==REGION_5PRIME ? "5prime" : "3prime");
    if (region != REGION_FULL) fprintf(stderr, "(%d bp)", region_len);
    fprintf(stderr, "\n");

    /* Allocate hash table */
    g_htable = calloc(HASH_SIZE, sizeof(Slot));
    if (!g_htable) { fprintf(stderr, "OOM\n"); return 1; }

    int n_kmers = build_kmer_hash(target_seq, target_len, kmer_size,
                                  region, region_len);
    if (n_kmers < 0) return 1;
    fprintf(stderr, "[%s/%s] k=%d  index: %d unique k-mers (fwd+RC)\n",
            sample, lane, kmer_size, n_kmers);

    uint64_t mask = (kmer_size < 32)
                  ? (((uint64_t)1 << (2*kmer_size)) - 1) : UINT64_MAX;

    long r2_sz = get_file_size(r2_path);

    gzFile fh1 = gzopen(r1_path, "rb");
    gzFile fh2 = gzopen(r2_path, "rb");
    if (!fh1 || !fh2) {
        fprintf(stderr, "[%s/%s] ERROR: cannot open input files\n", sample, lane);
        return 1;
    }
    gzbuffer(fh1, GZBUF_BYTES);
    gzbuffer(fh2, GZBUF_BYTES);

    FILE *out = fopen(outpath, "w");
    if (!out) {
        fprintf(stderr, "[%s/%s] ERROR: cannot open output: %s\n", sample, lane, outpath);
        return 1;
    }
    fprintf(out, "sample,lane,barcode,umi,n_kmer_hits,target_pos,strand,read_name,r2_seq\n");

    char h1[BUF_SIZE],s1[BUF_SIZE],p1[BUF_SIZE],q1[BUF_SIZE];
    char h2[BUF_SIZE],s2[BUF_SIZE],p2[BUF_SIZE],q2[BUF_SIZE];

    long n_reads=0, n_hits=0;
    time_t t_start=time(NULL), t_last=t_start;

    fprintf(stderr, "[%s/%s] Starting  R2 compressed: %.2f GB\n",
            sample, lane, r2_sz>0 ? (double)r2_sz/(1<<30) : 0.0);

    while (read_fastq_record(fh1,h1,s1,p1,q1) &&
           read_fastq_record(fh2,h2,s2,p2,q2))
    {
        n_reads++;
        int s2len = (int)strlen(s2);
        int first_pos, first_strand;
        int hits = scan_read(s2, s2len, kmer_size, mask, &first_pos, &first_strand);

        if (hits >= min_kmers) {
            n_hits++;
            int s1len    = (int)strlen(s1);
            int bc_end   = (barcode_len < s1len) ? barcode_len : s1len;
            int umi_start = bc_end;
            int umi_end   = umi_start + umi_len;
            if (umi_end > s1len) umi_end = s1len;

            const char *rn = (h2[0]=='@') ? h2+1 : h2;

            fprintf(out, "%s,%s,", sample, lane);
            fwrite(s1,            1, bc_end,             out); fputc(',',out);
            fwrite(s1+umi_start,  1, umi_end-umi_start,  out); fputc(',',out);
            fprintf(out, "%d,%d,%s,",
                    hits, first_pos, first_strand>0 ? "+" : "-");
            for (const char *c=rn; *c && *c!=' ' && *c!='\t'; c++) fputc(*c,out);
            fputc(',',out);
            csv_field(out, s2, s2len);
            fputc('\n',out);
        }

        if (n_reads % report_interval == 0) {
            time_t t_now = time(NULL);
            double elapsed  = difftime(t_now, t_start);
            double interval = difftime(t_now, t_last);
            double rps = (interval>0) ? (double)report_interval/interval : 0.0;
            t_last = t_now;

            double pct = 0.0;
            if (r2_sz > 0) {
                z_off_t pos = gztell(fh2);
                if (pos >= 0) pct = 100.0*(double)pos/(double)r2_sz;
            }
            char eta_buf[32];
            double eta = (pct>0.5) ? elapsed*(100.0/pct-1.0) : -1.0;
            fmt_time(eta, eta_buf, sizeof(eta_buf));

            fprintf(stderr,
                "[%s/%s]  %6.1fM reads | %5.1f%% | %5.0fk rd/s | ETA %s | hits: %ld (%.5f%%)\n",
                sample, lane,
                n_reads/1e6, pct, rps/1e3,
                eta_buf, n_hits, n_reads>0 ? 100.0*n_hits/n_reads : 0.0);
        }
    }

    double elapsed_total = difftime(time(NULL), t_start);
    char elapsed_buf[32];
    fmt_time(elapsed_total, elapsed_buf, sizeof(elapsed_buf));

    fprintf(stderr,
        "[%s/%s]  DONE | %.2fM reads | hits: %ld (%.5f%%) | Time: %s\n",
        sample, lane, n_reads/1e6, n_hits,
        n_reads>0 ? 100.0*n_hits/n_reads : 0.0, elapsed_buf);

    fclose(out);
    gzclose(fh1);
    gzclose(fh2);
    free(g_htable);
    return 0;
}
