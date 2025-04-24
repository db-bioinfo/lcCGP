import pandas as pd
import numpy as np
import json
from pathlib import Path
import argparse
from datetime import datetime

# Clinical thresholds as constants
MSI_THRESHOLD = 20.0  # Standard clinical threshold for MSI-H
MIN_COVERAGE_THRESHOLD = 100  # Minimum coverage for PASS
OPTIMAL_COVERAGE_THRESHOLD = 250  # Optimal coverage threshold
HIGH_QUALITY_SITES_THRESHOLD = 0.95  # 95% of sites should meet minimum coverage

def analyze_msi_data(sample_id, msi_file, dist_file, all_file, unstable_file):
    """
    Analyze MSI data from multiple input files and create structured output
    """
    # Read main MSI results
    with open(msi_file, 'r') as f:
        lines = f.readlines()
        total_sites, unstable_sites, percentage = lines[1].strip().split('\t')
        
    # Use standardized MSI threshold
    msi_status = "MSI-H" if float(percentage) >= MSI_THRESHOLD else "MSS"
    
    # Read distribution file for pattern analysis
    dist_df = pd.read_csv(dist_file, sep='\t')
    
    # Read all sites and unstable sites
    all_sites_df = pd.read_csv(all_file, sep='\t')
    unstable_df = pd.read_csv(unstable_file, sep='\t')
    
    # Analyze repeat patterns
    repeat_analysis = analyze_repeat_patterns(all_sites_df, unstable_df)
    
    # Analyze coverage
    coverage_analysis = analyze_coverage(all_sites_df, unstable_df)
    
    # Calculate quality status
    quality_status = determine_quality_status(coverage_analysis)
    
    # Create structured output
    analysis_result = {
        "sample_id": sample_id,
        "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "msi_status": {
            "classification": msi_status,
            "total_sites": int(total_sites),
            "unstable_sites": int(unstable_sites),
            "instability_percentage": float(percentage),
            "threshold_used": MSI_THRESHOLD
        },
        "quality_metrics": coverage_analysis,
        "quality_status": quality_status,
        "repeat_patterns": repeat_analysis,
        "notable_sites": get_notable_sites(unstable_df),
        "technical_details": {
            "tool": "MSIsensor-pro",
            "coverage_thresholds": {
                "minimum": MIN_COVERAGE_THRESHOLD,
                "optimal": OPTIMAL_COVERAGE_THRESHOLD
            },
            "instability_threshold": MSI_THRESHOLD/100  # Convert to decimal
        }
    }
    
    return analysis_result

def determine_quality_status(coverage_analysis):
    """
    Determine overall quality status based on coverage metrics
    """
    all_sites = coverage_analysis["all_sites"]
    status = {
        "overall": "PASS",
        "flags": [],
        "warnings": []
    }
    
    # Check mean coverage
    if all_sites["mean_coverage"] < MIN_COVERAGE_THRESHOLD:
        status["overall"] = "FAIL"
        status["flags"].append(f"Mean coverage ({all_sites['mean_coverage']:.1f}x) below minimum threshold ({MIN_COVERAGE_THRESHOLD}x)")
    elif all_sites["mean_coverage"] < OPTIMAL_COVERAGE_THRESHOLD:
        status["warnings"].append(f"Mean coverage ({all_sites['mean_coverage']:.1f}x) below optimal threshold ({OPTIMAL_COVERAGE_THRESHOLD}x)")
    
    # Check percentage of high-quality sites
    sites_fraction = all_sites["sites_above_min_coverage"] / all_sites["total_sites"]
    if sites_fraction < HIGH_QUALITY_SITES_THRESHOLD:
        status["overall"] = "FAIL"
        status["flags"].append(f"Only {sites_fraction*100:.1f}% of sites meet minimum coverage threshold (target: {HIGH_QUALITY_SITES_THRESHOLD*100}%)")
    
    return status

def analyze_coverage(all_df, unstable_df):
    """
    Analyze coverage statistics with updated thresholds
    """
    total_sites = len(all_df)
    return {
        "all_sites": {
            "mean_coverage": float(all_df['CovReads'].mean()),
            "median_coverage": float(all_df['CovReads'].median()),
            "min_coverage": int(all_df['CovReads'].min()),
            "max_coverage": int(all_df['CovReads'].max()),
            "sites_above_min_coverage": int((all_df['CovReads'] >= MIN_COVERAGE_THRESHOLD).sum()),
            "sites_above_optimal": int((all_df['CovReads'] >= OPTIMAL_COVERAGE_THRESHOLD).sum()),
            "total_sites": total_sites,
            "coverage_distribution": {
                f">{MIN_COVERAGE_THRESHOLD}x": float((all_df['CovReads'] >= MIN_COVERAGE_THRESHOLD).sum() / total_sites),
                f">{OPTIMAL_COVERAGE_THRESHOLD}x": float((all_df['CovReads'] >= OPTIMAL_COVERAGE_THRESHOLD).sum() / total_sites)
            }
        },
        "unstable_sites": {
            "mean_coverage": float(unstable_df['CovReads'].mean()),
            "median_coverage": float(unstable_df['CovReads'].median()),
            "min_coverage": int(unstable_df['CovReads'].min()),
            "max_coverage": int(unstable_df['CovReads'].max()),
            "sites_below_min_coverage": int((unstable_df['CovReads'] < MIN_COVERAGE_THRESHOLD).sum())
        }
    }

def analyze_repeat_patterns(all_df, unstable_df):
    """
    Analyze patterns in microsatellite repeats
    """
    # Analyze repeat unit types
    all_repeat_units = all_df['repeat_unit_bases'].value_counts()
    unstable_repeat_units = unstable_df['repeat_unit_bases'].value_counts()
    
    # Analyze repeat lengths
    all_df['repeat_length'] = all_df['repeat_times'].astype(int)
    unstable_df['repeat_length'] = unstable_df['repeat_times'].astype(int)
    
    # Calculate instability rates by repeat type
    repeat_type_stats = {}
    for repeat_unit in all_repeat_units.index:
        total_sites = all_repeat_units[repeat_unit]
        unstable_sites = unstable_repeat_units.get(repeat_unit, 0)
        instability_rate = (unstable_sites / total_sites) * 100 if total_sites > 0 else 0
        
        repeat_type_stats[repeat_unit] = {
            "total_sites": int(total_sites),
            "unstable_sites": int(unstable_sites),
            "instability_rate": float(instability_rate)
        }
    
    return {
        "repeat_unit_distribution": {
            "all_sites": all_repeat_units.to_dict(),
            "unstable_sites": unstable_repeat_units.to_dict()
        },
        "repeat_type_stats": repeat_type_stats,
        "repeat_length_stats": {
            "all_sites": {
                "mean": float(all_df['repeat_length'].mean()),
                "median": float(all_df['repeat_length'].median()),
                "std": float(all_df['repeat_length'].std())
            },
            "unstable_sites": {
                "mean": float(unstable_df['repeat_length'].mean()),
                "median": float(unstable_df['repeat_length'].median()),
                "std": float(unstable_df['repeat_length'].std())
            }
        }
    }

def get_notable_sites(unstable_df):
    """
    Identify and analyze notable unstable sites with additional metrics
    """
    # Sort by instability (pro_p value)
    top_sites = unstable_df.nlargest(5, 'pro_p')
    
    notable_sites = []
    for _, site in top_sites.iterrows():
        notable_sites.append({
            "location": f"{site['chromosome']}:{site['location']}",
            "repeat_unit": site['repeat_unit_bases'],
            "repeat_times": int(site['repeat_times']),
            "instability": float(site['pro_p']),
            "coverage": int(site['CovReads']),
            "quality_status": "PASS" if site['CovReads'] >= MIN_COVERAGE_THRESHOLD else "LOW_COVERAGE"
        })
    
    return notable_sites

def main():
    parser = argparse.ArgumentParser(description='Analyze MSI results and create structured output')
    parser.add_argument('--sample', required=True, help='Sample ID')
    parser.add_argument('--output', required=True, help='Output JSON file path')
    parser.add_argument('--msi', required=True, help='Main MSI results file')
    parser.add_argument('--dist', required=True, help='MSI distribution file')
    parser.add_argument('--all', required=True, help='All sites file')
    parser.add_argument('--unstable', required=True, help='Unstable sites file')
    
    args = parser.parse_args()
    
    try:
        # Run analysis
        result = analyze_msi_data(
            args.sample,
            args.msi,
            args.dist,
            args.all,
            args.unstable
        )
        
        # Save results
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        
        print(f"Analysis complete. Results saved to {args.output}")
        
        # Print quality status summary
        print("\nQuality Status Summary:")
        print(f"Overall Status: {result['quality_status']['overall']}")
        if result['quality_status']['flags']:
            print("\nQuality Flags:")
            for flag in result['quality_status']['flags']:
                print(f"- {flag}")
        if result['quality_status']['warnings']:
            print("\nWarnings:")
            for warning in result['quality_status']['warnings']:
                print(f"- {warning}")
            
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        raise

if __name__ == "__main__":
    main()
