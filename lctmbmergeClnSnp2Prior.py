#!/usr/bin/env python3
"""
Enhanced Variant Annotation Merger

This script merges annotations from a clnsnp file into a mutect file, 
matching variants by chromosome, position, and alleles with enhanced 
capabilities for handling complex variants and different representations.

Usage:
  python merge_variant_annotations.py <mutect_file> <clnsnp_file> <output_file>

Author: Clinical Bioinformatics Team
"""

import sys
import csv
import re
import logging
from collections import defaultdict

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# Gene name aliases have been removed as requested

def left_normalize_indel(ref, alt, pos):
    """
    Normalize indel representation by left-trimming common prefixes
    and right-trimming common suffixes.
    
    Args:
        ref (str): Reference allele
        alt (str): Alternate allele
        pos (int): Position
        
    Returns:
        tuple: (normalized_ref, normalized_alt, normalized_pos)
    """
    # Handle special cases
    if ref == alt:
        return ref, alt, pos
    
    if ref == '-':
        return '', alt, pos
    
    if alt == '-':
        return ref, '', pos
    
    # Left-trim common prefix
    original_length = len(ref)
    i = 0
    while i < min(len(ref), len(alt)) and ref[i] == alt[i]:
        i += 1
    
    if i > 0:
        ref = ref[i:]
        alt = alt[i:]
        pos += i
    
    # Right-trim common suffix
    j = 0
    while j < min(len(ref), len(alt)) and ref[-(j+1)] == alt[-(j+1)]:
        j += 1
    
    if j > 0:
        ref = ref[:-j]
        alt = alt[:-j]
    
    # If we've trimmed everything, represent as empty strings
    if not ref:
        ref = ''
    if not alt:
        alt = ''
    
    return ref, alt, pos

def parse_rs_id(rs_id):
    """Extract numeric RS ID from string with 'rs' prefix."""
    if not rs_id or rs_id == '.':
        return None
    match = re.search(r'rs(\d+)', rs_id)
    return match.group(1) if match else None

def normalize_chrom(chrom):
    """Normalize chromosome notation to ensure consistent matching."""
    if chrom.startswith('chr'):
        return chrom
    return f"chr{chrom}"

def get_gene_aliases(gene):
    """Get all possible aliases for a gene."""
    # Simply return the gene itself as we've removed the aliases functionality
    return [gene]

def find_best_match(mutect_variant, clnsnp_variants, gene_index, multi_allelic_variants):
    """
    Find the best matching variant in clnsnp for the given mutect variant.
    Uses multiple matching strategies with descending priority.
    
    Args:
        mutect_variant (dict): Variant from mutect file
        clnsnp_variants (list): List of all variants from clnsnp file
        gene_index (dict): Index of clnsnp variants by gene
        multi_allelic_variants (list): Expanded list of multi-allelic variants
    
    Returns:
        dict: Best matching variant from clnsnp, or None if no match found
    """
    chrom = mutect_variant['Chr']
    try:
        pos = int(mutect_variant['Start'])
    except ValueError:
        logger.warning(f"Invalid position in mutect variant: {mutect_variant['Start']}")
        pos = 0  # Default to 0 if position is invalid
        
    ref = mutect_variant['Ref']
    alt = mutect_variant['Alt']
    gene = mutect_variant['Ref.Gene']
    rs_id = parse_rs_id(mutect_variant.get('avsnp151', None))
    
    # If ref or alt is '-', convert to empty string for normalization
    if ref == '-':
        ref = ''
    if alt == '-':
        alt = ''
    
    # Get normalized representation
    norm_ref, norm_alt, norm_pos = left_normalize_indel(ref, alt, pos)
    
    # Get all possible gene names (including aliases)
    gene_aliases = get_gene_aliases(gene)
    
    # Strategy 1: Direct position and allele match (original logic)
    for v in clnsnp_variants:
        if (v['CHROM'] == chrom and 
            int(v['POS']) == pos and 
            v['REF'] == ref and 
            v['ALT'] == alt):
            return v
    
    # Strategy 2: Match normalized variants
    for v in clnsnp_variants:
        v_ref = v['REF']
        v_alt = v['ALT']
        
        # Skip multi-allelic variants here (they're handled separately)
        if ',' in v_alt:
            continue
            
        if v_ref == '-':
            v_ref = ''
        if v_alt == '-':
            v_alt = ''
            
        v_pos = int(v['POS'])
        v_norm_ref, v_norm_alt, v_norm_pos = left_normalize_indel(v_ref, v_alt, v_pos)
        
        if (v['CHROM'] == chrom and 
            v_norm_pos == norm_pos and 
            v_norm_ref == norm_ref and 
            v_norm_alt == norm_alt):
            return v
    
    # Strategy 3: Check multi-allelic variants
    for v in multi_allelic_variants:
        v_ref = v['REF']
        v_alt = v['ALT']
        
        if v_ref == '-':
            v_ref = ''
        if v_alt == '-':
            v_alt = ''
            
        v_pos = int(v['POS'])
        v_norm_ref, v_norm_alt, v_norm_pos = left_normalize_indel(v_ref, v_alt, v_pos)
        
        if (v['CHROM'] == chrom and 
            v_norm_pos == norm_pos and 
            v_norm_ref == norm_ref and 
            v_norm_alt == norm_alt):
            return v
    
    # Strategy 4: Match by RS ID if available (original logic)
    if rs_id:
        for v in clnsnp_variants:
            if v.get('RS', '').strip() == rs_id:
                return v
    
    # Strategy 5: Match by position with tolerance of ±1 (original logic)
    position_matches = []
    for offset in [-1, 0, 1]:
        for v in clnsnp_variants:
            if (v['CHROM'] == chrom and 
                int(v['POS']) == pos + offset and 
                v['ANN[0].GENE'] == gene):
                position_matches.append((abs(offset), v))
    
    if position_matches:
        # Sort by position offset (smallest first)
        position_matches.sort(key=lambda x: x[0])
        return position_matches[0][1]
    
    # Strategy 6: Match by position with larger tolerance for indels (±3)
    if len(ref) > 1 or len(alt) > 1 or ref == '-' or alt == '-':
        for offset in range(-3, 4):  # -3 to +3
            if offset in [-1, 0, 1]:  # Already checked these
                continue
                
            for v in clnsnp_variants:
                v_ref = v['REF']
                v_alt = v['ALT']
                
                # Skip multi-allelic variants
                if ',' in v_alt:
                    continue
                    
                # For indels, check if the variant is similar
                if (v['CHROM'] == chrom and 
                    abs(int(v['POS']) - pos) == abs(offset) and
                    (len(v_ref) > 1 or len(v_alt) > 1 or v_ref == '-' or v_alt == '-')):
                    # Check for gene match
                    if v['ANN[0].GENE'] == gene:
                        return v
    
    # Strategy 7: Match by gene in close proximity (±5 bp)
    # Note: Gene aliases functionality has been removed
    if gene in gene_index:
        candidates = gene_index[gene]
        for v in candidates:
            if v['CHROM'] == chrom and abs(int(v['POS']) - pos) <= 5:
                return v
    
    # Strategy 8: Match by HGVS notation if available (original logic)
    hgvs_c = None
    otherinfo = mutect_variant.get('Otherinfo', '')
    hgvs_match = re.search(r'c\.([A-Za-z0-9_>.+-]+)', otherinfo)
    if hgvs_match:
        hgvs_c = f"c.{hgvs_match.group(1)}"
        for v in clnsnp_variants:
            if v['ANN[0].GENE'] == gene and v['ANN[0].HGVS_C'] == hgvs_c:
                return v
    
    # No match found
    return None

def expand_multi_allelic_variants(clnsnp_variants):
    """
    Expand multi-allelic variants into individual variants for matching.
    
    Args:
        clnsnp_variants (list): List of variants from clnsnp file
        
    Returns:
        list: Expanded list of variants with one alternate allele each
    """
    expanded_variants = []
    
    for variant in clnsnp_variants:
        if ',' in variant['ALT']:
            # This is a multi-allelic variant
            alt_alleles = variant['ALT'].split(',')
            
            # Check if VAF is also comma-separated
            vaf_values = []
            if 'VAF' in variant and ',' in variant['VAF']:
                vaf_values = variant['VAF'].split(',')
            
            # Create a separate variant for each alternate allele
            for i, alt in enumerate(alt_alleles):
                new_variant = variant.copy()
                new_variant['ALT'] = alt
                
                # Update VAF if available
                if vaf_values and i < len(vaf_values):
                    new_variant['VAF'] = vaf_values[i]
                
                expanded_variants.append(new_variant)
        else:
            # Single-allelic variant, no expansion needed
            expanded_variants.append(variant)
    
    return expanded_variants

def main(mutect_file, clnsnp_file, output_file):
    # Columns to add from clnsnp
    columns_to_add = [
        "ANN[0].FEATUREID", 
        "ANN[0].EFFECT", 
        "ANN[0].HGVS_C", 
        "ANN[0].HGVS_P", 
        "Filter", 
        "CLNHGVS", 
        "VAF", 
        "RS"
    ]
    
    logger.info(f"Reading clnsnp file: {clnsnp_file}")
    
    # Read clnsnp file
    clnsnp_variants = []
    gene_index = defaultdict(list)  # Index clnsnp variants by gene
    
    with open(clnsnp_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            # Clean up potential Windows line endings in column names
            clean_row = {k.strip('\r'): v.strip('\r') for k, v in row.items()}
            clnsnp_variants.append(clean_row)
            gene = clean_row.get('ANN[0].GENE')
            if gene:
                gene_index[gene].append(clean_row)
    
    logger.info(f"Loaded {len(clnsnp_variants)} variants from clnsnp file")
    
    # Expand multi-allelic variants for matching
    multi_allelic_variants = expand_multi_allelic_variants(clnsnp_variants)
    logger.info(f"Expanded to {len(multi_allelic_variants)} variants after handling multi-allelic sites")
    
    # Read mutect file and create output file
    matched_count = 0
    unmatched_count = 0
    
    logger.info(f"Processing mutect file: {mutect_file}")
    logger.info(f"Writing output to: {output_file}")
    
    with open(mutect_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        
        # Create new header with added columns
        fieldnames = reader.fieldnames + columns_to_add
        
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        # Process each mutect variant
        for mutect_variant in reader:
            # Find best match in clnsnp
            best_match = find_best_match(
                mutect_variant, 
                clnsnp_variants, 
                gene_index,
                multi_allelic_variants
            )
            
            # Add columns from clnsnp match or use '.' for missing data
            for column in columns_to_add:
                if best_match and column in best_match:
                    # Clean up any Windows line endings
                    mutect_variant[column] = best_match[column].strip('\r')
                else:
                    mutect_variant[column] = "."
            
            # Track matching statistics
            if best_match:
                matched_count += 1
            else:
                unmatched_count += 1
                # Log unmatched variants for future improvement
                logger.debug(f"Unmatched: {mutect_variant['Chr']}:{mutect_variant['Start']} "
                           f"{mutect_variant['Ref']}>{mutect_variant['Alt']} ({mutect_variant['Ref.Gene']})")
            
            # Write the merged variant to output
            writer.writerow(mutect_variant)
    
    # Report statistics
    logger.info(f"Matching complete: {matched_count} matched, {unmatched_count} unmatched")
    logger.info(f"Match rate: {matched_count/(matched_count+unmatched_count)*100:.1f}%")
    logger.info(f"Merged annotations written to {output_file}")

def run_from_command_line():
    """Parse command line arguments and run the script."""
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <mutect_file> <clnsnp_file> <output_file>")
        sys.exit(1)
    
    mutect_file = sys.argv[1]
    clnsnp_file = sys.argv[2]
    output_file = sys.argv[3]
    
    main(mutect_file, clnsnp_file, output_file)

if __name__ == "__main__":
    run_from_command_line()
