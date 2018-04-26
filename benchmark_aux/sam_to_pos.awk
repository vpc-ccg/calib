function abs(v) {
    return v < 0 ? -v : v
}

BEGIN {
    F_QNAME = 1; # String [!-?A-~]{1,254} Query template NAME
    F_FLAG = 2; # Int [0,216-1] bitwise FLAG
    F_RNAME = 3; # String \*|[!-()+-<>-~][!-~]* Reference sequence NAME
    F_POS = 4; # Int [0,231-1] 1-based leftmost mapping POSition
    F_MAPQ = 5; # Int [0,28-1] MAPping Quality
    F_CIGAR = 6; # String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
    F_RNEXT = 7; # String \*|=|[!-()+-<>-~][!-~]* Ref. name of the mate/next read
    F_PNEXT = 8; # Int [0,231-1] Position of the mate/next read
    F_TLEN = 9; # Int [-231+1,231-1] observed Template LENgth
    F_SEQ = 10; # String \*|[A-Za-z=.]+ segment SEQuence
    F_QUAL = 11; # String [!-~]+ ASCII of Phred-scaled

    FLAG_MULTIPLE_SEQMENTS = 0x1; # template having multiple segments in sequencing
    FLAG_SEGMENTS_PROPERLY_ALIGNED = 0x2; # each segment properly aligned according to the aligner
    FLAG_SEGMENT_UNMAPPED = 0x4; # segment unmapped
    FLAG_NEXT_SEGMENT_UNMAPPED = 0x8; # next segment in the template unmapped
    FLAG_SEQ_REV_CMPL = 0x10; # SEQ being reverse complemented
    FLAG_NEXT_SEQ_REV_CMPL = 0x20; # SEQ of the next segment in the template being reverse complemented
    FLAG_FIRST_SEGMENT = 0x40; # the first segment in the template
    FLAG_LAST_SEGMENT = 0x80; # the last segment in the template
    FLAG_SECONDARY_MAPPING = 0x100; # secondary alignment
    FLAG_NOT_PASSING = 0x200; # not passing filters, such as platform/vendor quality controls
    FLAG_PCR_DUP = 0x400; # PCR or optical duplicate
    FLAG_SUPP_MAPPING = 0x800; # supplementary alignment
}

{
    if (substr($0,1,1) != "@") {
        # print and($F_FLAG, FLAG_SEGMENT_UNMAPPED)"\t"and($F_FLAG, FLAG_SEQ_REV_CMPL)
        printf $F_QNAME"\t";
        if ($F_POS == 0 || and($F_FLAG, FLAG_SEGMENT_UNMAPPED)) {
            printf unmapped_count--"\t"unmapped_count--"\t";
        } else if (and($F_FLAG, FLAG_SEQ_REV_CMPL)) {
            if (match($F_CIGAR, /[0-9]+S$/) < 0) {
                suffix_softclip_len = 0;
            } else {
                match($F_CIGAR, /[0-9]+S$/);
                suffix_softclip_len = int(substr($F_CIGAR, RSTART, RLENGTH-1));
            }
            printf $F_RNAME"\t"$F_PNEXT-$F_TLEN+suffix_softclip_len"\t";
        } else {
            if (match($F_CIGAR, /^S[0-9]+/) < 0){
                prefix_softclip_len = 0;
            } else {
                match($F_CIGAR, /^[0-9]+S/);
                prefix_softclip_len = int(substr($F_CIGAR, RSTART, RLENGTH-1));
            }
            printf $F_RNAME"\t"$F_POS-prefix_softclip_len"\t";
        }
        printf $F_SEQ"\n";#$0
    }
}
