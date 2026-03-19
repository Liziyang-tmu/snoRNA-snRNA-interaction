import sys
import argparse

def classify_snorna_with_overlap(input_file, output_base, threshold=10, min_overlap_ratio=0.5):

    counts = {
        'Total': 0, 'Mature': 0, 'Precursor': 0,
        'LowOverlap_Skipped': 0
    }
    
    mature_file = f"{output_base}_mature.bed"
    precursor_file = f"{output_base}_precursor.bed"
    
    with open(input_file, 'r') as f_in, \
         open(mature_file, 'w') as f_mature, \
         open(precursor_file, 'w') as f_precursor:
        

        header = "#read_chrom\tread_start\tread_end\tread_name\tscore\tread_strand\t" \
                "sno_chrom\tsno_start\tsno_end\tsno_name\tsno_score\tsno_strand\t" \
                "overlap_len\tread_len\tsno_len\tread_overlap_ratio\tleft_overhang\tright_overhang\tclassification\tnotes\n"
        f_mature.write(header)
        f_precursor.write(header)
        
        for line_num, line in enumerate(f_in, 1):
            if line.startswith('#') or line.strip() == '':
                continue
            
            cols = line.strip().split('\t')
            counts['Total'] += 1
            
            try:

                read_start = int(cols[1])
                read_end = int(cols[2])
                read_len = read_end - read_start  
                

                sno_start = int(cols[7])
                sno_end = int(cols[8])
                sno_len = sno_end - sno_start
                sno_strand = cols[11] 
                

                overlap_len = int(cols[-1])
                
            except (IndexError, ValueError) as e:
                print(f"Warning line {line_num}: Format error {e}")
                continue
            

            ratio_on_read = overlap_len / read_len if read_len > 0 else 0
            

            left_overhang = sno_start - read_start
            right_overhang = read_end - sno_end
            
            notes = []
            

            if ratio_on_read < min_overlap_ratio:
                counts['LowOverlap_Skipped'] += 1
                continue  
            

            is_precursor = (left_overhang > threshold or right_overhang > threshold)
            

            ratio_on_sno = overlap_len / sno_len if sno_len > 0 else 0
            
            if is_precursor:
                category = "Precursor"
                counts['Precursor'] += 1
                notes.append(f"boundary_extend(left={left_overhang},right={right_overhang})")
            else:
                category = "Mature"
                counts['Mature'] += 1
                notes.append(f"within_boundaries_or_short_overhang")
            
            notes.append(f"read_cov={ratio_on_read:.2f}")
            notes.append(f"sno_cov={ratio_on_sno:.2f}")
            
            notes_str = ";".join(notes)
            

            output_line = f"{cols[0]}\t{read_start}\t{read_end}\t{cols[3]}\t{cols[4]}\t{cols[5]}\t" \
                         f"{cols[6]}\t{sno_start}\t{sno_end}\t{cols[9]}\t{cols[10]}\t{cols[11]}\t" \
                         f"{overlap_len}\t{read_len}\t{sno_len}\t{ratio_on_read:.3f}\t" \
                         f"{left_overhang}\t{right_overhang}\t{category}\t{notes_str}\n"
            
            if category == "Mature":
                f_mature.write(output_line)
            else:
                f_precursor.write(output_line)
    
    print("\n" + "="*60)
    print("CLASSIFICATION SUMMARY (FIXED LOGIC)")
    print("="*60)
    print(f"Total reads processed: {counts['Total']:,}")
    print(f"Skipped (Low overlap): {counts['LowOverlap_Skipped']:,} (Read overlap ratio < {min_overlap_ratio})")
    print(f"Classified Mature:     {counts['Mature']:,}")
    print(f"Classified Precursor:  {counts['Precursor']:,}")
    print(f"\nDefinition used:")
    print(f"  - Precursor: Read extends > {threshold}bp outside snoRNA genomic boundaries")
    print(f"  - Mature:    Read is inside snoRNA (or extension <= {threshold}bp)")
    print("="*60)


def main():
    parser = argparse.ArgumentParser(description="Classify reads for snoRNA and crossing snoRNA 5' and 3' ends")
    parser.add_argument("-i", "--input", required=True, help="Input BED file (intersect -wo)")
    parser.add_argument("-o", "--output", required=True, help="Output base name")
    parser.add_argument("-t", "--threshold", type=int, default=10, help="Extension threshold bp (default: 10)")
    parser.add_argument("-r", "--overlap-ratio", type=float, default=0.5, 
                       help="Min ratio of READ length covered by snoRNA (default: 0.5)")
    
    args = parser.parse_args()
    classify_snorna_with_overlap(args.input, args.output, args.threshold, args.overlap_ratio)

if __name__ == "__main__":
    main()