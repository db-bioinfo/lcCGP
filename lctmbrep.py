#!/usr/bin/env python3
import sys
import os
import json
import csv
import re
from datetime import datetime

def read_tsv_file(file_path):
    """Read TSV file and return headers and data rows"""
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)  # Get the first row as headers
        rows = []
        for i, row in enumerate(reader):
            if row and len(row) >= len(headers):  # Skip empty or invalid rows
                # Create a dictionary for each row
                variant = {headers[j]: row[j] if j < len(row) else '' for j in range(len(headers))}
                # Add rank based on position
                variant['rank'] = i + 1
                rows.append(variant)
    return headers, rows

def parse_tmb_report(file_path):
    """Parse TMB report file and extract relevant information"""
    tmb_data = {
        'tmb_score': None,
        'effective_genome_size': None,
        'min_depth': None,
        'total_variants': None,
        'variants_after_filters': None
    }
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
            # Extract TMB score
            tmb_match = re.search(r'TMB=\s*(\d+\.\d+)', content)
            if tmb_match:
                tmb_data['tmb_score'] = float(tmb_match.group(1))
                
            # Extract Effective Genome Size
            genome_match = re.search(r'Effective Genome Size=\s*(\d+)', content)
            if genome_match:
                tmb_data['effective_genome_size'] = int(genome_match.group(1))
                
            # Extract minimum depth
            depth_match = re.search(r'minDepth=\s*(\d+)', content)
            if depth_match:
                tmb_data['min_depth'] = int(depth_match.group(1))
                
            # Extract total variants
            total_variants_match = re.search(r'Total number of variants=\s*(\d+)', content)
            if total_variants_match:
                tmb_data['total_variants'] = int(total_variants_match.group(1))
                
            # Extract variants after filters
            filtered_variants_match = re.search(r'Variants after filters=\s*(\d+)', content)
            if filtered_variants_match:
                tmb_data['variants_after_filters'] = int(filtered_variants_match.group(1))
    
    except Exception as e:
        print(f"Error parsing TMB report: {e}")
    
    return tmb_data

def parse_tumor_purity(file_path):
    """Parse tumor purity CSV file and extract relevant information"""
    purity_data = {
        'purity': None,
        'ploidy': None
    }
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                # Just get the first row's data
                purity_data['purity'] = float(row.get('Purity', 0)) if row.get('Purity') else None
                purity_data['ploidy'] = float(row.get('Ploidy', 0)) if row.get('Ploidy') else None
                break  # We only need the first row
    
    except Exception as e:
        print(f"Error parsing tumor purity file: {e}")
    
    return purity_data

def parse_stats_file(file_path):
    """Parse Stats.txt file and extract average coverage"""
    stats_data = {
        'average_coverage': None
    }
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
            
            # Extract Average Coverage
            coverage_match = re.search(r'AVERAGE_COVERAGE:(\d+\.\d+)', content)
            if coverage_match:
                stats_data['average_coverage'] = float(coverage_match.group(1))
    
    except Exception as e:
        print(f"Error parsing Stats file: {e}")
    
    return stats_data

def generate_html(variants_file, tmb_report_file=None, tumor_purity_file=None, stats_file=None, output_path=None):
    """Generate HTML report from TSV file with optional TMB and tumor purity data"""
    # Determine output filename if not specified
    if not output_path:
        output_path = os.path.splitext(variants_file)[0] + "_report.html"
    
    # Get the sample name from the file path (for IGV linking)
    full_filename = os.path.basename(os.path.splitext(variants_file)[0])
    # Extract just the sample ID (numbers before any underscore)
    sample_id_match = re.match(r'^(\d+)', full_filename)
    sample_name = sample_id_match.group(1) if sample_id_match else full_filename
    
    # Read the variants TSV file
    try:
        headers, variants = read_tsv_file(variants_file)
    except Exception as e:
        print(f"Error reading variants TSV file: {e}")
        return False
    
    # Parse TMB report if provided
    tmb_data = {}
    if tmb_report_file and os.path.exists(tmb_report_file):
        tmb_data = parse_tmb_report(tmb_report_file)
    
    # Parse tumor purity if provided
    purity_data = {}
    if tumor_purity_file and os.path.exists(tumor_purity_file):
        purity_data = parse_tumor_purity(tumor_purity_file)
    
    # Parse stats file if provided
    stats_data = {}
    if stats_file and os.path.exists(stats_file):
        stats_data = parse_stats_file(stats_file)
    
    # Convert data to JSON for embedding in HTML
    variants_json = json.dumps(variants)
    tmb_json = json.dumps(tmb_data)
    purity_json = json.dumps(purity_data)
    stats_json = json.dumps(stats_data)
    
    # Current date for the report
    current_date = datetime.now().strftime("%B %d, %Y")
    
    # Generate HTML content with embedded data
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Genetic Variant Analysis Report</title>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css" rel="stylesheet">
    <style>
        /* Color Palette as specified */
        :root {{
            /* Primary Colors (Main Branding) */
            --deep-blue: #06274b;
            --teal: #2CA6A4;
            --soft-gray: #F4F4F4;
            
            /* Accent Colors (For Highlights, Buttons, and Graphs) */
            --orange: #F39237;
            --green: #4CAF50;
            --red: #E74C3C;
            --purple: #9C27B0;
            --lime: #a6ce4d; /* New lime green color for icons and bars */
            
            /* Typography & Background */
            --dark-text: #212121;
            --white-bg: #FFFFFF;
            
            /* Additional UI Colors */
            --light-border: #E0E0E0;
            --medium-gray: #757575;
            --light-teal: rgba(44, 166, 164, 0.1);
            
            /* Custom colors for badges */
            --pathogenic-red: #E74C3C;
            --likely-pathogenic-red: #F07470; /* Lighter red for likely pathogenic */
            --conflicting-red: #FF8C00; /* dark orange for conflicting */
            
            /* ACMG code colors */
            --path-color: #E74C3C;
            --benign-color: #4CAF50;
            
            /* TMB Classification Colors */
            --tmb-low: #4CAF50;      /* Green */
            --tmb-int-low: #FFC107;  /* Yellow */
            --tmb-int-high: #FF9800; /* Orange */
            --tmb-high: #E74C3C;     /* Red */
        }}
        
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
        }}
        
        body {{
            background-color: var(--soft-gray);
            color: var(--dark-text);
            line-height: 1.6;
        }}
        
        .container {{
            max-width: 100%;
            margin: 0 auto;
            padding: 20px;
        }}
        
        /* Header Styles - REDUCED BY 15% to match WES report */
        .report-header {{
            background-color: var(--deep-blue);
            color: var(--white-bg);
            padding: 8.5px 12px; /* Reduced from 10px */
            border-radius: 10px 10px 0 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: flex;
            justify-content: center;
            align-items: center;
            width: 23%;
            max-width: 231px;
            margin-left: 0;
            margin-right: auto;
            height: 45px; /* Reduced from ~51px (15% reduction) */
        }}
        
        .logo img {{
            height: 29px; /* Reduced from 34px (15% reduction) */
            width: auto;
        }}
        
        /* Info Card */
        .info-card {{
            background-color: var(--light-teal);
            border-left: 4px solid var(--teal);
            padding: 10px 15px;
            margin-bottom: 20px;
            border-radius: 4px;
        }}
        
        .info-card h3 {{
            color: var(--deep-blue);
            margin-bottom: 8px;
            font-size: 16px;
        }}
        
        .info-details {{
            display: flex;
            flex-direction: column;
            gap: 6px;
        }}
        
        .info-item {{
            display: flex;
            align-items: flex-start;
            gap: 6px;
        }}
        
        .info-icon {{
            color: var(--teal);
            width: 16px;
            text-align: center;
            font-size: 12px;
            margin-top: 3px;
        }}
        
        .info-text {{
            font-size: 13px;
            color: var(--dark-text);
            line-height: 1.4;
        }}
        
        /* NEW: TMB Metrics Section - Styled like WES Coverage Metrics */
        .coverage-section {{
            background-color: var(--white-bg);
            padding: 10px; /* Reduced from 20px */
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 15px; /* Reduced from 20px */
        }}
        
        .coverage-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 8px; /* Reduced from 15px */
            border-bottom: 1px solid var(--light-border);
            padding-bottom: 5px; /* Reduced from 10px */
        }}
        
        .coverage-header h3 {{
            color: var(--deep-blue);
            font-size: 14px; /* Reduced from 16px */
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        
        .coverage-metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 8px; /* Reduced from 15px */
        }}
        
        .metrics-card {{
            background-color: var(--soft-gray);
            border-radius: 8px;
            padding: 8px; /* Reduced from 15px */
            box-shadow: 0 1px 2px rgba(0,0,0,0.05);
            transition: transform 0.2s, box-shadow 0.2s;
            height: 55px; /* Reduced height - was approx 80-90px */
        }}
        
        .metrics-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 3px 6px rgba(0,0,0,0.1);
        }}
        
        .metrics-title {{
            font-size: 11px; /* Reduced from 12px */
            color: var(--medium-gray);
            margin-bottom: 3px; /* Reduced from 5px */
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        
        .metrics-title i {{
            color: var(--lime); /* Changed from teal to lime green */
        }}
        
        .metrics-value {{
            font-size: 15px; /* Reduced from 18px */
            font-weight: 600;
            color: var(--deep-blue);
        }}
        
        /* Controls */
        .controls {{
            background-color: var(--white-bg);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            display: flex;
            flex-direction: column;
        }}
        
        .search-controls {{
            display: flex;
            gap: 15px;
            align-items: center;
            margin-bottom: 15px;
        }}
        
        .search-box {{
            flex-grow: 1;
            padding: 10px 15px;
            border: 1px solid var(--light-border);
            border-radius: 4px;
            font-size: 14px;
        }}
        
        .refresh-button {{
            display: inline-flex;
            align-items: center;
            gap: 5px;
            background-color: var(--teal);
            color: white;
            border: none;
            padding: 6px 12px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background-color 0.2s;
            white-space: nowrap;
        }}
        
        .refresh-button:hover {{
            background-color: #249391;
        }}
        
        /* Stats and Chart Container */
        .stats-chart-container {{
            background-color: var(--white-bg);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
            display: grid;
            grid-template-columns: 3fr 2fr;
            gap: 20px;
        }}
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
        }}
        
        .metric-card {{
            background-color: var(--soft-gray);
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        
        .metric-value {{
            font-size: 1.8rem;
            font-weight: bold;
            color: var(--deep-blue);
        }}
        
        .metric-label {{
            color: var(--medium-gray);
            font-size: 0.9rem;
            margin-top: 5px;
        }}
        
        .chart-container {{
            height: 100%;
            min-height: 250px;
            display: flex;
            justify-content: center;
            align-items: center;
        }}
        
        /* Variants Table */
        .variants-container {{
            background-color: var(--white-bg);
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            overflow: auto; /* Enable horizontal scrolling */
            margin-bottom: 20px;
            position: relative;
        }}
        
        .variants-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
            min-width: 1600px; /* Adjusted for merged columns */
            table-layout: fixed;
        }}
        
        .variants-table th,
        .variants-table td {{
            padding: 12px 15px;
            text-align: left;
            position: relative;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        
        /* Tooltip styles */
        .variants-table td.has-tooltip:hover {{
            overflow: visible !important;
            white-space: normal !important;
            z-index: 20;
        }}
        
        /* UPDATED: Improved Table Headers with report26 style */
        .variants-table th {{
            background: linear-gradient(180deg, var(--deep-blue) 0%, #2C5282 100%);
            color: var(--white-bg);
            font-weight: 600;
            cursor: pointer;
            position: sticky;
            top: 0;
            z-index: 10;
            transition: background 0.2s ease;
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        
        .variants-table th::after {{
            content: '';
            position: absolute;
            right: 10px;
            top: 50%;
            transform: translateY(-50%);
            opacity: 0.5;
            transition: opacity 0.3s ease;
        }}
        
        .variants-table th[data-direction="asc"]::after {{
            content: '\\f0d8';
            font-family: 'Font Awesome 6 Free';
            font-weight: 900;
            opacity: 1;
        }}
        
        .variants-table th[data-direction="desc"]::after {{
            content: '\\f0d7';
            font-family: 'Font Awesome 6 Free';
            font-weight: 900;
            opacity: 1;
        }}
        
        .variants-table th:hover {{
            background: linear-gradient(180deg, #2C5282 0%, #1A365D 100%);
        }}
        
        .variants-table tbody tr:nth-child(even) {{
            background-color: var(--soft-gray);
        }}
        
        .variants-table tbody tr:hover {{
            background-color: rgba(44, 166, 164, 0.05);
        }}
        
        .variants-table td {{
            border-bottom: 1px solid var(--light-border);
            overflow: hidden;
        }}

        /* Table column widths */
        .col-rank {{ width: 55px; }}
        .col-variant {{ width: 80px; }}
        .col-rs {{ width: 35px; }} /* RS column moved after variant */
        .col-gene {{ width: 100px; }}
        .col-effect {{ width: 145px; }}
        .col-clnhgvs {{ width: 120px; }}
        .col-hgvs {{ width: 140px; }} /* Merged HGVS column with increased width */
        .col-clinvar {{ width: 215px; }}
        .col-acmg {{ width: 165px; }}
        .col-acmg-codes {{ width: 130px; }}
        .col-af-gnomad {{ width: 105px; }}
        .col-af {{ width: 80px; }} /* Replaces zygosity */
        .col-phenotype {{ width: 60px; }}
        .col-disease {{ width: 60px; }}
        .col-database {{ width: 140px; }} /* Reduced width for icon-based database column */
        .col-igv {{ width: 50px; }} /* IGV column */
        .col-filter {{ width: 80px; }} /* Filter column */
        
        /* Special handling for variant code */
        .variant-code {{
            font-family: 'Courier New', monospace;
            background-color: rgba(30, 58, 95, 0.1);
            padding: 2px 5px;
            border-radius: 3px;
            font-size: 13px;
            max-width: 160px;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
            display: inline-block;
        }}
        
        /* Rank Styling */
        .rank-bubble {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 32px;
            height: 32px;
            background-color: var(--deep-blue);
            color: var(--white-bg);
            border-radius: 50%;
            font-weight: bold;
        }}
        
        /* Gene Styling */
        .gene-link {{
            color: var(--teal);
            text-decoration: none;
            font-weight: 600;
        }}
        
        .gene-link:hover {{
            text-decoration: underline;
        }}
        
        /* UPDATED: Badge Styling for ACMG and Clinical Significance with report26 style */
        .badge {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            padding: 4px 10px;
            border-radius: 50px;
            font-size: 12px;
            font-weight: 600;
            color: var(--white-bg);
            text-align: center;
            transition: all 0.2s ease;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            white-space: nowrap;
            text-transform: none;
            letter-spacing: 0;
        }}
        
        .badge:hover {{
            transform: translateY(-1px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.15);
        }}
        
        .badge-pathogenic {{
            background: linear-gradient(90deg, var(--pathogenic-red) 0%, #DC2626 100%);
        }}
        
        .badge-likely-pathogenic {{
            background: linear-gradient(90deg, var(--likely-pathogenic-red) 0%, #EF4444 100%);
        }}
        
        .badge-uncertain {{
            background: linear-gradient(90deg, var(--teal) 0%, #0891B2 100%);
        }}
        
        .badge-benign, .badge-likely-benign {{
            background: linear-gradient(90deg, var(--green) 0%, #059669 100%);
        }}
        
        .badge-conflicting {{
            background: linear-gradient(90deg, var(--orange) 0%, #D97706 100%);
        }}
        
        .badge-other {{
            background: linear-gradient(90deg, var(--purple) 0%, #7C3AED 100%);
        }}
        
        /* UPDATED: Database links styling with report26 style */
        .database-links {{
            display: flex;
            gap: 6px;
            flex-wrap: wrap;
            justify-content: flex-start;
        }}
        
        /* Icon-style database links */
        .db-icon-link {{
            display: inline-flex;
            align-items: center;
            justify-content: center;
            width: 28px;
            height: 28px;
            border-radius: 50%;
            font-size: 12px;
            font-weight: bold;
            text-decoration: none;
            transition: all 0.2s ease;
            margin: 2px;
            position: relative;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }}
        
        /* Database-specific icon styles */
        .db-icon-link.clinvar {{
            background-color: #E2F1FF;
            color: #0366D6;
        }}
        
        .db-icon-link.varsome {{
            background-color: #E4F8E9;
            color: #28A745;
        }}
        
        .db-icon-link.franklin {{
            background-color: #5087a3;
            color: #ffffff;
        }}
        
        .db-icon-link.omim {{
            background-color: #7e49a3;
            color: #ffffff;
        }}
        
        .db-icon-link.orphanet {{
            background-color: #94a3c2;
            color: #ffffff;
        }}
        
        /* IGV link styling */
        .db-icon-link.igv {{
            background-color: #d4d0a3;
            color: #ffffff;
        }}
        
        .db-icon-link:hover {{
            transform: translateY(-3px) scale(1.1);
            box-shadow: 0 4px 8px rgba(0,0,0,0.15);
        }}
        
        .db-icon-link.clinvar:hover {{ background-color: #0366D6; color: white; }}
        .db-icon-link.varsome:hover {{ background-color: #28A745; color: white; }}
        .db-icon-link.franklin:hover {{ background-color: #F39237; color: white; }}
        .db-icon-link.omim:hover {{ background-color: #6518B4; color: white; }}
        .db-icon-link.orphanet:hover {{ background-color: #D81B60; color: white; }}
        .db-icon-link.igv:hover {{ background-color: #3A0CA3; color: white; }}
        
        /* Pagination */
        .pagination {{
            display: flex;
            justify-content: center;
            margin-top: 20px;
            gap: 5px;
        }}
        
        .pagination button {{
            background-color: var(--white-bg);
            border: 1px solid var(--light-border);
            color: var(--deep-blue);
            padding: 8px 12px;
            cursor: pointer;
            border-radius: 4px;
            transition: all 0.2s;
        }}
        
        .pagination button:hover {{
            background-color: var(--light-teal);
        }}
        
        .pagination button.active {{
            background-color: var(--teal);
            color: var(--white-bg);
            border-color: var(--teal);
        }}
        
        /* Footer */
        .footer {{
            background-color: var(--white-bg);
            padding: 20px;
            border-radius: 0 0 10px 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-top: 20px;
            text-align: center;
            color: var(--medium-gray);
            font-size: 0.9rem;
        }}
        
        /* Not available text */
        .not-available {{
            color: var(--medium-gray);
            font-style: italic;
        }}
        
        /* ACMG Criteria - NEW STYLING OPTION 3 */
        .acmg-container {{
            display: flex;
            flex-wrap: wrap;
            gap: 8px;
        }}
        
        .acmg-code {{
            display: inline-block;
            padding: 2px 6px;
            font-family: monospace;
            font-size: 12px;
            background-color: var(--soft-gray);
            border-radius: 3px;
        }}
        
        .acmg-path::before {{
            content: '▼ ';
            color: var(--path-color);
        }}
        
        .acmg-benign::before {{
            content: '▲ ';
            color: var(--benign-color);
        }}
        
        /* Info icon styling for the data cells */
        .data-info-icon {{
            font-size: 16px;
            color: var(--lime); /* Changed from teal to lime green */
            cursor: pointer;
            text-align: center;
            transition: transform 0.2s;
        }}
        
        .data-info-icon:hover {{
            transform: scale(1.2);
        }}
        
        /* Center info icon in the cell */
        .center-icon {{
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%;
        }}
        
        /* Center N/A text in cells */
        .center-na {{
            display: flex;
            justify-content: center;
            align-items: center;
            height: 100%;
            color: var(--medium-gray);
            font-style: italic;
        }}
        
        /* HGVS stacked styles */
        .hgvs-container {{
            display: flex;
            flex-direction: column;
            gap: 5px;
        }}
        
        .hgvs-item {{
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
        
        .hgvs-item .label {{
            font-weight: bold;
            margin-right: 4px;
            font-size: 12px;
            color: var(--deep-blue);
        }}
        
        /* Fixed tooltip styles - IMPROVED FOR WRAPPING LONG CONTENT */
        .fixed-tooltip {{
            display: none;
            position: absolute;
            background-color: #1A8582; 
            color: white;
            padding: 8px 12px;
            border-radius: 4px;
            font-size: 15px;
            font-weight: bold;
            z-index: 1000;
            max-width: 350px; /* Increased from 300px */
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
            pointer-events: none;
            white-space: normal; /* Ensure wrapping */
            line-height: 1.4;
            overflow-wrap: break-word; /* Force words to break */
            word-wrap: break-word;
            word-break: break-word; /* Additional word breaking for long strings */
        }}
        
        /* Copy notification tooltip */
        .copy-tooltip {{
            position: absolute;
            background-color: #333;
            color: white;
            padding: 5px 10px;
            border-radius: 4px;
            font-size: 12px;
            z-index: 1100;
            box-shadow: 0 2px 8px rgba(0,0,0,0.2);
            opacity: 0;
            transition: opacity 0.3s ease;
            pointer-events: none;
        }}
        
        /* Responsive Design */
        @media (max-width: 1200px) {{
            .container {{
                padding: 10px;
            }}
            
            .stats-chart-container {{
                grid-template-columns: 1fr;
            }}
            
            .chart-container {{
                min-height: 200px;
            }}
            
            .coverage-metrics-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
            
            .coverage-bars {{
                grid-template-columns: repeat(2, 1fr);
            }}
        }}
        
        @media (max-width: 768px) {{
            .report-header {{
                width: 41%;
                max-width: none;
            }}
            
            .metrics-grid {{
                grid-template-columns: 1fr 1fr;
            }}
            
            .coverage-metrics-grid {{
                grid-template-columns: 1fr;
            }}
            
            .coverage-bars {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <!-- Fixed tooltip container that follows cursor -->
    <div id="fixedTooltip" class="fixed-tooltip"></div>
    
    <!-- Copy notification tooltip -->
    <div id="copyTooltip" class="copy-tooltip">Copied!</div>
    
    <div class="container">
        <!-- Report Header - MODIFIED to match WES report -->
        <div class="report-header">
            <div class="logo">
                <img src="logo.png" alt="Company Logo" style="height: 29px; width: auto;">
            </div>
        </div>
        
        <!-- Info Card - Updated with more professional formatting -->
        <div class="info-card">
            <h3>Sample Information</h3>
            <div class="info-details">
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-dna"></i></div>
                    <div class="info-text"><strong>TMB Results</strong> • Sample ID: {sample_name}</div>
                </div>
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-calendar-alt"></i></div>
                    <div class="info-text"><strong>Analysis Date:</strong> {current_date}</div>
                </div>
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-clipboard-list"></i></div>
                    <div class="info-text">This report contains prioritized genetic variants with TMB (Tumor Mutational Burden) analysis. Variants are classified according to ACMG guidelines and analyzed for clinical relevance.</div>
                </div>
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-sort-amount-down"></i></div>
                    <div class="info-text">Variants are sorted by clinical significance, with potentially pathogenic variants appearing first.</div>
                </div>
            </div>
        </div>
        
        <!-- TMB Analysis Section - Styled like WES Coverage Metrics -->
        <div class="coverage-section" id="tmbDashboard">
            <div class="coverage-header">
                <h3><i class="fas fa-dna"></i> Tumor Mutational Burden</h3>
            </div>
            
            <div class="coverage-metrics-grid">
                <!-- TMB Score -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-calculator"></i> TMB Score</div>
                    <div class="metrics-value" id="tmbScoreValue">-- mut/Mb</div>
                </div>
                
                <!-- Total Variants -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-file-alt"></i> Total Variants</div>
                    <div class="metrics-value" id="totalVariants">--</div>
                </div>
                
                <!-- Variants After Filters -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-filter"></i> Variants After Filters</div>
                    <div class="metrics-value" id="variantsAfterFilters">--</div>
                </div>
                
                <!-- Average Coverage -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-layer-group"></i> Average Coverage</div>
                    <div class="metrics-value" id="averageCoverage">-- X</div>
                </div>
                
                <!-- Tumor Purity -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-vial"></i> Tumor Purity</div>
                    <div class="metrics-value" id="tumorPurity">--</div>
                </div>
                
                <!-- Tumor Ploidy -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-dna"></i> Tumor Ploidy</div>
                    <div class="metrics-value" id="tumorPloidy">--</div>
                </div>
                
                <!-- Effective Genome Size -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-ruler"></i> Effective Genome</div>
                    <div class="metrics-value" id="effectiveGenomeSize">--</div>
                </div>
                
                <!-- Min Sequencing Depth -->
                <div class="metrics-card">
                    <div class="metrics-title"><i class="fas fa-microscope"></i> Min Depth</div>
                    <div class="metrics-value" id="minDepth">--</div>
                </div>
            </div>
        </div>
        
        <!-- Search and Reset View Button - Modified to match WES report style -->
        <div class="controls">
            <div class="search-controls">
                <button id="refreshButton" class="refresh-button"><i class="fas fa-sync-alt"></i> Reset View</button>
                <input type="text" id="searchInput" class="search-box" placeholder="Search for genes, variants, or phenotypes...">
            </div>
        </div>
        
        <!-- Variants Table -->
        <div class="variants-container">
            <table class="variants-table" id="variantsTable">
                <thead>
                    <tr>
                        <th class="col-rank" data-sort="rank">RANK</th>
                        <th class="col-variant" data-sort="variant">Variant</th>
                        <th class="col-rs" data-sort="rs">rs</th>
                        <th class="col-gene" data-sort="gene">Gene</th>
                        <th class="col-effect" data-sort="effect">Effect</th>
                        <th class="col-clnhgvs" data-sort="clnhgvs">CLNHGVS</th>
                        <th class="col-hgvs" data-sort="hgvs">HGVS</th>
                        <th class="col-clinvar" data-sort="clinvar">Clinvar</th>
                        <th class="col-acmg" data-sort="acmg">ACMG</th>
                        <th class="col-acmg-codes" data-sort="acmgCodes">ACMG Codes</th>
                        <th class="col-af-gnomad" data-sort="afGnomad">AF GnomAD</th>
                        <th class="col-af" data-sort="af">AF</th>
                        <th class="col-phenotype" data-sort="phenotype">P <i class="fas fa-notes-medical"></i></th>
                        <th class="col-disease" data-sort="disease">D <i class="fas fa-notes-medical"></i></th>
                        <th class="col-database" data-sort="database">Database</th>
                        <th class="col-igv" data-sort="igv">IGV</th>
                        <th class="col-filter" data-sort="filter">Filter</th>
                    </tr>
                </thead>
                <tbody id="variantsTableBody">
                    <!-- Variants will be added here dynamically -->
                </tbody>
            </table>
        </div>
        
        <!-- Pagination -->
        <div class="pagination" id="pagination">
            <!-- Pagination will be dynamically generated -->
        </div>
        
        <!-- Footer -->
        <div class="footer">
            <p>Generated with Genomic Variant Analysis Pipeline</p>
            <p><small>For research and clinical use. Sort columns by clicking headers.</small></p>
        </div>
    </div>

    <script>
        // Store the sample name for IGV links
        const sampleName = "{sample_name}";
        
        // Embedded variant data from TSV file
        const allVariants = {variants_json};
        
        // Embedded TMB data
        const tmbData = {tmb_json};
        
        // Embedded Tumor Purity data
        const purityData = {purity_json};
        
        // Embedded Stats data
        const statsData = {stats_json};
        
        document.addEventListener('DOMContentLoaded', function() {{
            // Initialize variables
            let filteredVariants = [...allVariants];
            let currentPage = 1;
            const rowsPerPage = 15;
            let originalVariants = [...allVariants]; // Keep original ordering for reset
            
            // Initialize TMB Dashboard if data is available
            initializeTMBDashboard();
            
            // Function to initialize TMB Dashboard with data
            function initializeTMBDashboard() {{
                // Check if TMB data is available
                if (Object.keys(tmbData).length > 0) {{
                    // Set TMB Score
                    const tmbScoreValue = document.getElementById('tmbScoreValue');
                    
                    if (tmbData.tmb_score !== null) {{
                        const score = tmbData.tmb_score;
                        tmbScoreValue.textContent = score.toFixed(2) + ' mut/Mb';
                    }}
                    
                    // Set Total Variants
                    const totalVariants = document.getElementById('totalVariants');
                    if (tmbData.total_variants !== null) {{
                        totalVariants.textContent = tmbData.total_variants.toLocaleString();
                    }}
                    
                    // Set Variants After Filters
                    const variantsAfterFilters = document.getElementById('variantsAfterFilters');
                    if (tmbData.variants_after_filters !== null) {{
                        variantsAfterFilters.textContent = tmbData.variants_after_filters.toLocaleString();
                    }}
                    
                    // Set Effective Genome Size
                    const effectiveGenomeSize = document.getElementById('effectiveGenomeSize');
                    if (tmbData.effective_genome_size !== null) {{
                        effectiveGenomeSize.textContent = tmbData.effective_genome_size.toLocaleString();
                    }}
                    
                    // Set Minimum Depth
                    const minDepth = document.getElementById('minDepth');
                    if (tmbData.min_depth !== null) {{
                        minDepth.textContent = tmbData.min_depth;
                    }}
                }}
                
                // Check if Stats data is available for Average Coverage
                if (Object.keys(statsData).length > 0) {{
                    // Set Average Coverage
                    const averageCoverage = document.getElementById('averageCoverage');
                    if (statsData.average_coverage !== null) {{
                        averageCoverage.textContent = statsData.average_coverage.toFixed(1) + ' X';
                    }}
                }}
                
                // Check if Tumor Purity data is available
                if (Object.keys(purityData).length > 0) {{
                    // Set Tumor Purity
                    const tumorPurity = document.getElementById('tumorPurity');
                    if (purityData.purity !== null) {{
                        tumorPurity.textContent = (purityData.purity * 100).toFixed(1) + '%';
                    }}
                    
                    // Set Tumor Ploidy
                    const tumorPloidy = document.getElementById('tumorPloidy');
                    if (purityData.ploidy !== null) {{
                        tumorPloidy.textContent = purityData.ploidy.toFixed(2);
                    }}
                }}
                
                // If no TMB data is available, hide the dashboard
                if (Object.keys(tmbData).length === 0 && Object.keys(purityData).length === 0 && Object.keys(statsData).length === 0) {{
                    const tmbDashboard = document.getElementById('tmbDashboard');
                    if (tmbDashboard) {{
                        tmbDashboard.style.display = 'none';
                    }}
                }}
            }}
            
            // Priority mappings for sorting
            const acmgPriority = {{
                'Pathogenic': 1,
                'Likely pathogenic': 2,
                'Uncertain significance': 3,
                'Likely benign': 4,
                'Benign': 5
            }};
            
            const clinvarPriority = {{
                'Pathogenic': 1,
                'Pathogenic/Likely_pathogenic': 2,
                'Likely_pathogenic': 3,
                'Conflicting Classification': 4,
                'Conflicting_classifications_of_pathogenicity': 4,
                'Uncertain_significance': 5,
                'Likely_benign': 6,
                'Benign': 7,
                'Benign/Likely_benign': 8,
                'drug_response': 9,
                'risk_factor': 10,
                'protective': 11,
                'association': 12,
                'Affects': 13,
                'other': 14,
                'UNK': 15
            }};
            
            const tierPriority = {{
                'Tier 1 - High Priority': 1,
                'Tier 2 - Medium Priority': 2
            }};
            
            // Function to copy text to clipboard
            function copyToClipboard(text) {{
                // Create a temporary textarea element
                const textarea = document.createElement('textarea');
                textarea.value = text;
                
                // Make the textarea invisible and add it to the document
                textarea.style.position = 'absolute';
                textarea.style.left = '-9999px';
                document.body.appendChild(textarea);
                
                // Select and copy the text
                textarea.select();
                document.execCommand('copy');
                
                // Remove the textarea
                document.body.removeChild(textarea);
                
                // Show copy feedback tooltip
                const copyTooltip = document.getElementById('copyTooltip');
                
                // Position the tooltip near the mouse cursor
                const mousePosition = {{
                    x: window.event ? window.event.clientX : 0,
                    y: window.event ? window.event.clientY : 0
                }};
                
                // Position it slightly above the cursor
                copyTooltip.style.left = (mousePosition.x + 10) + 'px';
                copyTooltip.style.top = (mousePosition.y - 30) + 'px';
                
                // Show the tooltip
                copyTooltip.style.opacity = '1';
                
                // Hide it after a short delay
                setTimeout(() => {{
                    copyTooltip.style.opacity = '0';
                }}, 1500);
                
                return true;
            }}
            
            // Function to add tooltip to a cell
            function addTooltipToCell(cell, content, makeClickable = false) {{
                if (!content || content === 'N/A' || content === '.') return;
                
                cell.setAttribute('data-full-content', content);
                
                // Show tooltip on mouseenter
                cell.addEventListener('mouseenter', function(e) {{
                    const tooltip = document.getElementById('fixedTooltip');
                    tooltip.textContent = this.getAttribute('data-full-content');
                    tooltip.style.display = 'block';
                }});
                
                // Hide tooltip on mouseleave
                cell.addEventListener('mouseleave', function() {{
                    document.getElementById('fixedTooltip').style.display = 'none';
                }});
                
                // If the cell should be clickable for copying content
                if (makeClickable) {{
                    cell.style.cursor = 'pointer';
                    cell.addEventListener('click', function() {{
                        const textToCopy = this.getAttribute('data-full-content');
                        if (textToCopy) {{
                            copyToClipboard(textToCopy);
                        }}
                    }});
                }}
            }}
            
            // Set up tooltip functionality
            document.addEventListener('mousemove', function(e) {{
                const tooltip = document.getElementById('fixedTooltip');
                if (tooltip.style.display === 'block') {{
                    // Find the cell that triggered the tooltip
                    let targetCell = document.elementFromPoint(e.clientX, e.clientY);
                    while (targetCell && targetCell.tagName !== 'TD') {{
                        targetCell = targetCell.parentElement;
                    }}
                    
                    // Position based on column class
                    if (targetCell && (
                        targetCell.classList.contains('col-phenotype') || 
                        targetCell.classList.contains('col-disease') || 
                        targetCell.classList.contains('col-database') ||
                        targetCell.classList.contains('col-igv')
                    )) {{
                        // Position tooltip to the left of the cursor
                        tooltip.style.left = (e.pageX - tooltip.offsetWidth - 15) + 'px';
                    }} else if (targetCell && targetCell.classList.contains('col-rs')) {{
                        // Position tooltip to the right of the cursor for RS column
                        tooltip.style.left = (e.pageX + 15) + 'px';
                    }} else {{
                        // Default: Position tooltip to the right of the cursor
                        tooltip.style.left = (e.pageX + 15) + 'px';
                    }}
                    
                    // Position tooltip above cursor if near bottom of screen
                    const windowHeight = window.innerHeight;
                    const tooltipHeight = tooltip.offsetHeight;
                    if (e.clientY + tooltipHeight + 25 > windowHeight) {{
                        tooltip.style.top = (e.pageY - tooltipHeight - 10) + 'px';
                    }} else {{
                        tooltip.style.top = (e.pageY + 15) + 'px';
                    }}
                }}
            }});
            
            // Functions for variant processing
            function getAcmgBadgeClass(classification) {{
                if (!classification) return '';
                
                classification = classification.toLowerCase();
                if (classification.includes('pathogenic') && !classification.includes('likely')) {{
                    return 'badge-pathogenic';
                }} else if (classification.includes('likely pathogenic')) {{
                    return 'badge-likely-pathogenic';
                }} else if (classification.includes('uncertain') || classification.includes('vus')) {{
                    return 'badge-uncertain';
                }} else if (classification.includes('likely benign')) {{
                    return 'badge-likely-benign';
                }} else if (classification.includes('benign')) {{
                    return 'badge-benign';
                }}
                
                return '';
            }}
            
            function getClinicvarBadgeClass(classification) {{
                if (!classification || classification === 'UNK') return '';
                
                classification = classification.toLowerCase();
                
                // 1. Check for conflicting classifications - should be red just like pathogenic
                if (classification.includes('conflicting')) {{
                    return 'badge-conflicting';
                }}
                // 2. Check for pathogenic (most severe) - red
                else if (classification.includes('pathogenic') && 
                        !classification.includes('likely_pathogenic') && 
                        !classification.includes('likely pathogenic')) {{
                    return 'badge-pathogenic';
                }} 
                // 3. Check for likely pathogenic - lighter red
                else if (classification.includes('likely_pathogenic') || 
                        classification.includes('likely pathogenic')) {{
                    return 'badge-likely-pathogenic';
                }} 
                // 4. Check for uncertain significance or VUS
                else if (classification.includes('uncertain') || classification.includes('vus')) {{
                    return 'badge-uncertain';
                }} 
                // 5. Check for likely benign
                else if (classification.includes('likely_benign') || 
                        classification.includes('likely benign')) {{
                    return 'badge-likely-benign';
                }} 
                // 6. Check for benign
                else if (classification.includes('benign')) {{
                    return 'badge-benign';
                }}
                // 7. Other classifications
                else {{
                    return 'badge-other';
                }}
            }}
            
            function formatVariant(chrom, pos, ref, alt) {{
                // Handle special case for very long sequences
                if (ref && ref.length > 20) {{
                    ref = ref.substring(0, 10) + '...' + ref.substring(ref.length - 10);
                }}
                if (alt && alt.length > 20) {{
                    alt = alt.substring(0, 10) + '...' + alt.substring(alt.length - 10);
                }}
                return `${{chrom}}:${{pos}}${{ref}}>${{alt}}`;
            }}
            
            function formatVariantWithTooltip(chrom, pos, ref, alt) {{
                // Create full representation for tooltip
                const fullVariant = `${{chrom}}:${{pos}}${{ref}}>${{alt}}`;
                
                // Create shortened version for display
                let displayVariant = fullVariant;
                if ((ref && ref.length > 20) || (alt && alt.length > 20)) {{
                    if (ref && ref.length > 20) {{
                        ref = ref.substring(0, 10) + '...' + ref.substring(ref.length - 10);
                    }}
                    if (alt && alt.length > 20) {{
                        alt = alt.substring(0, 10) + '...' + alt.substring(alt.length - 10);
                    }}
                    displayVariant = `${{chrom}}:${{pos}}${{ref}}>${{alt}}`;
                }}
                
                return {{ full: fullVariant, display: displayVariant }};
            }}
            
            // Updated function to get Effect from ANN[0].EFFECT first, then ExonicFunc.refGene
            // With modifications to replace * with space and remove 'variant'
            function getEffect(variant) {{
                let effectText = '';
                if (variant['ANN[0].EFFECT'] && variant['ANN[0].EFFECT'] !== '.') {{
                    effectText = variant['ANN[0].EFFECT'];
                }} else {{
                    effectText = variant['ExonicFunc.refGene'] && variant['ExonicFunc.refGene'] !== 'unknown' ? variant['ExonicFunc.refGene'] : 'N/A';
                }}
                
                // Replace * with space and remove 'variant'
                if (effectText && effectText !== 'N/A') {{
                    effectText = effectText.replace(/\*/g, ' ').replace(/_/g, ' ');
                    effectText = effectText.replace(/variant/g, '').trim();
                }}
                
                return effectText;
            }}
            
            function formatPhenotype(orphaText) {{
                if (!orphaText) return '';
                
                const phenotypes = [];
                const orphaEntries = orphaText.split('~');
                
                orphaEntries.forEach(entry => {{
                    const parts = entry.split('|');
                    if (parts.length > 1) {{
                        phenotypes.push(parts[1]);
                    }}
                }});
                
                return phenotypes.join(', ');
            }}
            
            // Extract Orphanet IDs from Orpha field
            function getOrphanetIDs(orphaText) {{
                if (!orphaText) return [];
                
                const orphaIds = [];
                const orphaEntries = orphaText.split('~');
                
                orphaEntries.forEach(entry => {{
                    const parts = entry.split('|');
                    if (parts.length > 0) {{
                        // Extract ID from format like "ORPHA:12345"
                        const idMatch = parts[0].match(/ORPHA:(\\d+)/);
                        if (idMatch && idMatch[1]) {{
                            orphaIds.push(idMatch[1]);
                        }}
                    }}
                }});
                
                return orphaIds;
            }}
            
            // Get disease information from CLNDN field
            function getDisease(variant) {{
                // Try different possible column names for CLNDN
                const possibleKeys = [
                    'CLNDN',
                    'clndn',
                    'Clndn'
                ];
                
                for (const key of possibleKeys) {{
                    if (variant[key] && variant[key] !== '.' && variant[key] !== 'N/A') {{
                        // Process CLNDN field - replace pipes and split by vertical bars if needed
                        return variant[key].replace(/\\|/g, ', ').replace(/;/g, ', ');
                    }}
                }}
                
                return '';
            }}
            
            // Function to get ClinVar value, trying different possible column names
            function getClinvarValue(variant) {{
                // Check various possible column names for ClinVar data
                const possibleKeys = [
                    'clinvar: Clinvar',
                    'clinvar:Clinvar',
                    'clinvar: Clinvar ',
                    'clinvar',
                    'Clinvar',
                    'CLINVAR',
                    'CLNSIG'
                ];
                
                for (const key of possibleKeys) {{
                    if (variant[key] && variant[key] !== '.' && variant[key] !== 'N/A') {{
                        // Remove "clinvar:" prefix if present
                        let value = variant[key];
                        if (value.startsWith('clinvar:')) {{
                            value = value.substring('clinvar:'.length).trim();
                        }}
                        if (value.startsWith('clinvar: ')) {{
                            value = value.substring('clinvar: '.length).trim();
                        }}
                        
                        // Keep original Uncertain_significance value
                        if (value === 'Uncertain_significance') {{
                            value = 'Uncertain_significance';
                        }}
                        // Replace long Conflicting classification string
                        else if (value.startsWith('Conflicting_classifications_of_pathogenicity')) {{
                            value = 'Conflicting Classification';
                        }}
                        
                        return value;
                    }}
                }}
                
                return null;
            }}
            
            // Function to generate database links for the variant
            function generateDatabaseLinks(variant) {{
                const links = {{
                    clinvar: '',
                    varsome: '',
                    franklin: '',
                    omim: '',       // Added for OMIM
                    orphanet: ''    // Added for Orphanet
                }};
                
                // Get RS ID for links
                const rsId = variant.avsnp151 && variant.avsnp151 !== '.' ? variant.avsnp151 : null;
                
                // Get HGVS for links - now using CLNHGVS
                const hgvs = variant.CLNHGVS && variant.CLNHGVS !== '.' ? variant.CLNHGVS : null;
                
                // Get gene for links
                const gene = variant['Ref.Gene'] || '';
                
                // Create ClinVar link - UPDATED to prioritize HGVS
                if (hgvs) {{
                    // Prioritize HGVS for ClinVar links
                    links.clinvar = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${{encodeURIComponent(hgvs)}}`;
                }} else if (rsId) {{
                    // Fall back to rsId if HGVS is not available
                    links.clinvar = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${{rsId}}`;
                }} else {{
                    // Create a link using chromosome, position, reference, and alternate alleles
                    const chr = variant.Chr || '';
                    const pos = variant.Start || '';
                    links.clinvar = `https://www.ncbi.nlm.nih.gov/clinvar/?term=${{chr}}:${{pos}}`;
                }}
                
                // Create VarSome link
                if (rsId) {{
                    links.varsome = `https://varsome.com/variant/hg19/${{rsId}}`;
                }} else if (hgvs) {{
                    links.varsome = `https://varsome.com/variant/hg19/${{encodeURIComponent(hgvs)}}`;
                }} else {{
                    // Format for VarSome using chr-pos-ref-alt format
                    const chr = variant.Chr || '';
                    const pos = variant.Start || '';
                    const ref = variant.Ref || '';
                    const alt = variant.Alt || '';
                    links.varsome = `https://varsome.com/variant/hg19/${{chr}}-${{pos}}-${{ref}}-${{alt}}`;
                }}
                
                // Create Franklin link - ALWAYS USE chr-pos-ref-alt FORMAT
                const chr = variant.Chr || '';
                const pos = variant.Start || '';
                const ref = variant.Ref || '';
                const alt = variant.Alt || '';
                links.franklin = `https://franklin.genoox.com/clinical-db/variant/snp/chr${{chr}}-${{pos}}-${{ref}}-${{alt}}`;
                
                // Create OMIM link
                if (variant.OMIM && variant.OMIM !== '.' && variant.OMIM !== '') {{
                    // Take the first OMIM number if there are multiple separated by semicolons
                    const omimId = variant.OMIM.split(';')[0].trim();
                    if (omimId) {{
                        links.omim = `https://www.omim.org/entry/${{omimId}}`;
                    }}
                }}
                
                // Create Orphanet link
                if (variant.OrphaNumber && variant.OrphaNumber !== '.' && variant.OrphaNumber !== '') {{
                    // Take the first Orphanet number if there are multiple separated by semicolons
                    const orphaId = variant.OrphaNumber.split(';')[0].trim();
                    if (orphaId) {{
                        links.orphanet = `https://www.orpha.net/en/disease/detail/${{orphaId}}?name=${{orphaId}}&mode=orpha`;
                    }}
                }}
                
                return links;
            }}
            
            // Function to generate IGV link for a variant
            function generateIGVLink(variant) {{
                // Use the exact IGV HTML file path format without location hash
                // The IGV.html file is created by the pipeline
                return `${{sampleName}}.IGV.html`;
            }}
            
            // Reset button functionality
            document.getElementById('refreshButton').addEventListener('click', function() {{
                // Reset the variants to original state
                filteredVariants = [...originalVariants];
                
                // Clear search input
                document.getElementById('searchInput').value = '';
                
                // Reset column headers (clear sorting indicators)
                const tableHeaders = document.querySelectorAll('#variantsTable th');
                tableHeaders.forEach(header => {{
                    header.removeAttribute('data-direction');
                }});
                
                // Reset to first page and update display
                currentPage = 1;
                renderTable(currentPage);
                renderPagination();
            }});
            
            // Render table with pagination
            function renderTable(page) {{
                const startIndex = (page - 1) * rowsPerPage;
                const endIndex = startIndex + rowsPerPage;
                const displayedVariants = filteredVariants.slice(startIndex, endIndex);
                
                const tbody = document.getElementById('variantsTableBody');
                tbody.innerHTML = '';
                
                // If no variants to display
                if (displayedVariants.length === 0) {{
                    tbody.innerHTML = `
                        <tr>
                            <td colspan="17" style="text-align: center; padding: 20px;">
                                No variants match your search criteria.
                            </td>
                        </tr>
                    `;
                    return;
                }}
                
                // Render each variant
                displayedVariants.forEach((variant, index) => {{
                    const row = document.createElement('tr');
                    
                    // RANK
                    const rankCell = document.createElement('td');
                    rankCell.className = 'col-rank';
                    const rankBubble = document.createElement('div');
                    rankBubble.className = 'rank-bubble';
                    rankBubble.textContent = startIndex + index + 1;
                    rankCell.appendChild(rankBubble);
                    // Add tooltip for rank
                    addTooltipToCell(rankCell, `Rank: ${{startIndex + index + 1}}`);
                    row.appendChild(rankCell);
                    
                    // MODIFIED: Variant - replace with info icon
                    const variantCell = document.createElement('td');
                    variantCell.className = 'col-variant';
                    
                    const formattedVariant = formatVariantWithTooltip(
                        variant.Chr, 
                        variant.Start, 
                        variant.Ref, 
                        variant.Alt
                    );
                    
                    // Create info icon container for centered alignment
                    const iconContainer = document.createElement('div');
                    iconContainer.className = 'center-icon';
                    
                    // Create info icon
                    const infoIcon = document.createElement('i');
                    infoIcon.className = 'fas fa-info-circle data-info-icon';
                    iconContainer.appendChild(infoIcon);
                    
                    // Add tooltip to the info icon with clickable option
                    addTooltipToCell(variantCell, formattedVariant.full, true);
                    
                    variantCell.appendChild(iconContainer);
                    row.appendChild(variantCell);
                    
                    // MODIFIED: RS - replace with info icon (moved here before gene)
                    const rsCell = document.createElement('td');
                    rsCell.className = 'col-rs';
                    const rsValue = variant.avsnp151 || 'N/A';
                    
                    if (rsValue === 'N/A' || rsValue === '.') {{
                        // Create centered N/A container
                        const naContainer = document.createElement('div');
                        naContainer.className = 'center-na';
                        naContainer.textContent = 'N/A';
                        rsCell.appendChild(naContainer);
                    }} else {{
                        // Create info icon container for centered alignment
                        const iconContainer = document.createElement('div');
                        iconContainer.className = 'center-icon';
                        
                        // Create info icon
                        const infoIcon = document.createElement('i');
                        infoIcon.className = 'fas fa-info-circle data-info-icon';
                        iconContainer.appendChild(infoIcon);
                        
                        // Add tooltip for rs with clickable option
                        addTooltipToCell(rsCell, rsValue, true);
                        
                        rsCell.appendChild(iconContainer);
                    }}
                    row.appendChild(rsCell);
                    
                    // Gene
                    const geneCell = document.createElement('td');
                    geneCell.className = 'col-gene';
                    const geneLink = document.createElement('a');
                    geneLink.href = `https://www.genecards.org/cgi-bin/carddisp.pl?gene=${{variant['Ref.Gene']}}`;
                    geneLink.target = '_blank';
                    geneLink.className = 'gene-link';
                    geneLink.textContent = variant['Ref.Gene'];
                    geneCell.appendChild(geneLink);
                    // Add tooltip for gene
                    addTooltipToCell(geneCell, variant['Ref.Gene']);
                    row.appendChild(geneCell);
                    
                    // Effect - always with tooltip
                    const effectCell = document.createElement('td');
                    effectCell.className = 'col-effect';
                    const effectText = getEffect(variant);
                    effectCell.textContent = effectText || 'N/A';
                    
                    if (!effectText || effectText === 'N/A') {{
                        effectCell.className += ' not-available';
                    }} else {{
                        addTooltipToCell(effectCell, effectText);
                    }}
                    row.appendChild(effectCell);
                    
                    // CLNHGVS - always with tooltip
                    const clnhgvsCell = document.createElement('td');
                    clnhgvsCell.className = 'col-clnhgvs';
                    const clnhgvsText = variant.CLNHGVS || 'N/A';
                    clnhgvsCell.textContent = clnhgvsText;
                    
                    if (clnhgvsText === 'N/A' || clnhgvsText === '.') {{
                        clnhgvsCell.className += ' not-available';
                    }} else {{
                        addTooltipToCell(clnhgvsCell, clnhgvsText);
                    }}
                    row.appendChild(clnhgvsCell);
                    
                    // NEW: Merged HGVS Column (combining HGVS C and HGVS P)
                    const hgvsCell = document.createElement('td');
                    hgvsCell.className = 'col-hgvs';
                    
                    const hgvsCText = variant['ANN[0].HGVS_C'] || 'N/A';
                    const hgvsPText = variant['ANN[0].HGVS_P'] || 'N/A';
                    
                    // Container for stacked HGVS values
                    const hgvsContainer = document.createElement('div');
                    hgvsContainer.className = 'hgvs-container';
                    
                    // HGVS C entry
                    const hgvsCItem = document.createElement('div');
                    hgvsCItem.className = 'hgvs-item';
                    
                    const hgvsCLabel = document.createElement('span');
                    hgvsCLabel.className = 'label';
                    hgvsCLabel.textContent = 'C:';
                    hgvsCItem.appendChild(hgvsCLabel);
                    
                    const hgvsCValue = document.createTextNode(hgvsCText !== 'N/A' ? hgvsCText : 'N/A');
                    hgvsCItem.appendChild(hgvsCValue);
                    
                    if (hgvsCText !== 'N/A' && hgvsCText !== '.') {{
                        addTooltipToCell(hgvsCItem, hgvsCText);
                    }} else {{
                        hgvsCItem.classList.add('not-available');
                    }}
                    
                    hgvsContainer.appendChild(hgvsCItem);
                    
                    // HGVS P entry
                    const hgvsPItem = document.createElement('div');
                    hgvsPItem.className = 'hgvs-item';
                    
                    const hgvsPLabel = document.createElement('span');
                    hgvsPLabel.className = 'label';
                    hgvsPLabel.textContent = 'P:';
                    hgvsPItem.appendChild(hgvsPLabel);
                    
                    const hgvsPValue = document.createTextNode(hgvsPText !== 'N/A' ? hgvsPText : 'N/A');
                    hgvsPItem.appendChild(hgvsPValue);
                    
                    if (hgvsPText !== 'N/A' && hgvsPText !== '.') {{
                        addTooltipToCell(hgvsPItem, hgvsPText);
                    }} else {{
                        hgvsPItem.classList.add('not-available');
                    }}
                    
                    hgvsContainer.appendChild(hgvsPItem);
                    hgvsCell.appendChild(hgvsContainer);
                    
                    row.appendChild(hgvsCell);
                    
                    // UPDATED: Clinvar with badge styling and updated naming
                    const clinvarCell = document.createElement('td');
                    clinvarCell.className = 'col-clinvar';
                    const clinvarVal = getClinvarValue(variant);
                    
                    if (clinvarVal && clinvarVal !== 'UNK') {{
                        const badge = document.createElement('span');
                        badge.className = `badge ${{getClinicvarBadgeClass(clinvarVal)}}`;
                        badge.textContent = clinvarVal;
                        clinvarCell.appendChild(badge);
                        // Tooltip for Clinvar removed as requested
                    }} else {{
                        clinvarCell.textContent = 'N/A';
                        clinvarCell.className += ' not-available';
                    }}
                    row.appendChild(clinvarCell);
                    
                    // UPDATED: ACMG with badge styling and VUS text update
                    const acmgCell = document.createElement('td');
                    acmgCell.className = 'col-acmg';
                    if (variant.ACMG && variant.ACMG !== '.') {{
                        const badge = document.createElement('span');
                        badge.className = `badge ${{getAcmgBadgeClass(variant.ACMG)}}`;
                        
                        // Keep original "Uncertain significance" text
                        let acmgText = variant.ACMG;
                        
                        badge.textContent = acmgText;
                        acmgCell.appendChild(badge);
                        // Tooltip for ACMG removed as requested
                    }} else {{
                        acmgCell.textContent = 'N/A';
                        acmgCell.className += ' not-available';
                    }}
                    row.appendChild(acmgCell);
                    
                    // ACMG Codes - NEW MONOCHROME WITH CHARACTER PREFIX STYLING
                    const acmgCodesCell = document.createElement('td');
                    acmgCodesCell.className = 'col-acmg-codes';
                    
                    if (variant['ACMG Criteria'] && variant['ACMG Criteria'] !== '.') {{
                        const criteriaContainer = document.createElement('div');
                        criteriaContainer.className = 'acmg-container';
                        
                        // Split the criteria codes
                        const codes = variant['ACMG Criteria'].split(',').map(code => code.trim());
                        
                        codes.forEach(code => {{
                            if (code) {{
                                const codeElement = document.createElement('span');
                                codeElement.className = 'acmg-code';
                                codeElement.textContent = code;
                                
                                // Add appropriate class based on code type
                                if (code.startsWith('P')) {{
                                    codeElement.className += ' acmg-path';
                                }} else if (code.startsWith('B')) {{
                                    codeElement.className += ' acmg-benign';
                                }}
                                
                                criteriaContainer.appendChild(codeElement);
                            }}
                        }});
                        
                        acmgCodesCell.appendChild(criteriaContainer);
                    }} else {{
                        // Create centered N/A container
                        const naContainer = document.createElement('div');
                        naContainer.className = 'center-na';
                        naContainer.textContent = 'N/A';
                        acmgCodesCell.appendChild(naContainer);
                    }}
                    
                    row.appendChild(acmgCodesCell);
                    
                    // AF GnomAD - no tooltip
                    const afGnomadCell = document.createElement('td');
                    afGnomadCell.className = 'col-af-gnomad';
                    if (variant['Freq_gnomAD_genome_ALL'] && variant['Freq_gnomAD_genome_ALL'] !== '.') {{
                        afGnomadCell.textContent = variant['Freq_gnomAD_genome_ALL'];
                        // Removed tooltip for AF GnomAD
                    }} else {{
                        afGnomadCell.textContent = 'N/A';
                        afGnomadCell.className += ' not-available';
                    }}
                    row.appendChild(afGnomadCell);
                    
                    // AF - replaces VAF/Zygosity - no tooltip
                    const afCell = document.createElement('td');
                    afCell.className = 'col-af';
                    const afText = variant.VAF || 'N/A';
                    afCell.textContent = afText;
                    
                    if (afText === 'N/A' || afText === '.') {{
                        afCell.className += ' not-available';
                    }}
                    row.appendChild(afCell);
                    
                    // MODIFIED: Phenotype - replace with info icon
                    const phenotypeCell = document.createElement('td');
                    phenotypeCell.className = 'col-phenotype';
                    const phenotypeText = formatPhenotype(variant.Orpha);
                    
                    if (!phenotypeText) {{
                        // Create centered N/A container
                        const naContainer = document.createElement('div');
                        naContainer.className = 'center-na';
                        naContainer.textContent = 'N/A';
                        phenotypeCell.appendChild(naContainer);
                    }} else {{
                        // Create info icon container for centered alignment
                        const iconContainer = document.createElement('div');
                        iconContainer.className = 'center-icon';
                        
                        // Create info icon
                        const infoIcon = document.createElement('i');
                        infoIcon.className = 'fas fa-info-circle data-info-icon';
                        iconContainer.appendChild(infoIcon);
                        
                        // Add tooltip for phenotype with clickable option
                        addTooltipToCell(phenotypeCell, phenotypeText, true);
                        
                        phenotypeCell.appendChild(iconContainer);
                    }}
                    row.appendChild(phenotypeCell);
                    
                    // MODIFIED: Disease - replace with info icon
                    const diseaseCell = document.createElement('td');
                    diseaseCell.className = 'col-disease';
                    const diseaseText = getDisease(variant);
                    
                    if (!diseaseText) {{
                        // Create centered N/A container
                        const naContainer = document.createElement('div');
                        naContainer.className = 'center-na';
                        naContainer.textContent = 'N/A';
                        diseaseCell.appendChild(naContainer);
                    }} else {{
                        // Create info icon container for centered alignment
                        const iconContainer = document.createElement('div');
                        iconContainer.className = 'center-icon';
                        
                        // Create info icon
                        const infoIcon = document.createElement('i');
                        infoIcon.className = 'fas fa-info-circle data-info-icon';
                        iconContainer.appendChild(infoIcon);
                        
                        // Add tooltip for disease with clickable option
                        addTooltipToCell(diseaseCell, diseaseText, true);
                        
                        diseaseCell.appendChild(iconContainer);
                    }}
                    row.appendChild(diseaseCell);
                    
                    // UPDATED: Database with Font Awesome icons like report26
                    const databaseCell = document.createElement('td');
                    databaseCell.className = 'col-database';
                    
                    // Generate database links
                    const dbLinks = generateDatabaseLinks(variant);
                    
                    // Create container for links
                    const linksContainer = document.createElement('div');
                    linksContainer.className = 'database-links';
                    
                    // Create ClinVar icon link with Font Awesome
                    const clinvarLink = document.createElement('a');
                    clinvarLink.href = dbLinks.clinvar;
                    clinvarLink.className = 'db-icon-link clinvar';
                    clinvarLink.innerHTML = '<i class="fas fa-flask"></i>';
                    clinvarLink.target = '_blank';
                    addTooltipToCell(clinvarLink, 'ClinVar Database');
                    linksContainer.appendChild(clinvarLink);
                    
                    // Create VarSome icon link with Font Awesome
                    const varsomeLink = document.createElement('a');
                    varsomeLink.href = dbLinks.varsome;
                    varsomeLink.className = 'db-icon-link varsome';
                    varsomeLink.innerHTML = '<i class="fas fa-dna"></i>';
                    varsomeLink.target = '_blank';
                    addTooltipToCell(varsomeLink, 'VarSome Database');
                    linksContainer.appendChild(varsomeLink);
                    
                    // Create Franklin icon link with Font Awesome
                    const franklinLink = document.createElement('a');
                    franklinLink.href = dbLinks.franklin;
                    franklinLink.className = 'db-icon-link franklin';
                    franklinLink.innerHTML = '<i class="fas fa-vial"></i>';
                    franklinLink.target = '_blank';
                    addTooltipToCell(franklinLink, 'Franklin Database');
                    linksContainer.appendChild(franklinLink);
                    
                    // Create OMIM icon link if available with Font Awesome
                    if (dbLinks.omim) {{
                        const omimLink = document.createElement('a');
                        omimLink.href = dbLinks.omim;
                        omimLink.className = 'db-icon-link omim';
                        omimLink.innerHTML = '<i class="fas fa-book-medical"></i>';
                        omimLink.target = '_blank';
                        addTooltipToCell(omimLink, 'OMIM Database');
                        linksContainer.appendChild(omimLink);
                    }}
                    
                    // Create Orphanet icon link if available with Font Awesome
                    if (dbLinks.orphanet) {{
                        const orphanetLink = document.createElement('a');
                        orphanetLink.href = dbLinks.orphanet;
                        orphanetLink.className = 'db-icon-link orphanet';
                        orphanetLink.innerHTML = '<i class="fas fa-notes-medical"></i>';
                        orphanetLink.target = '_blank';
                        addTooltipToCell(orphanetLink, 'Orphanet Database');
                        linksContainer.appendChild(orphanetLink);
                    }}
                    
                    databaseCell.appendChild(linksContainer);
                    row.appendChild(databaseCell);
                    
                    // IGV cell
                    const igvCell = document.createElement('td');
                    igvCell.className = 'col-igv';
                    
                    // Only create IGV link if variant rank is <= 300
                    // The IGV report only contains the first 300 variants
                    if (variant.rank <= 300) {{
                        // Create container for IGV link
                        const igvContainer = document.createElement('div');
                        igvContainer.className = 'center-icon';
                        
                        // Create IGV link with Font Awesome
                        const igvLink = document.createElement('a');
                        igvLink.href = generateIGVLink(variant);
                        igvLink.className = 'db-icon-link igv';
                        igvLink.style.height = '28px'; // Force same height as other database icons
                        igvLink.style.width = '28px';  // Force same width as other database icons
                        igvLink.innerHTML = '<i class="fas fa-eye"></i>';
                        igvLink.target = '_blank';
                        addTooltipToCell(igvLink, 'View in IGV Browser');
                        
                        igvContainer.appendChild(igvLink);
                        igvCell.appendChild(igvContainer);
                    }} else {{
                        // Create centered N/A container
                        const naContainer = document.createElement('div');
                        naContainer.className = 'center-na';
                        naContainer.textContent = 'N/A';
                        igvCell.appendChild(naContainer);
                    }}
                    
                    row.appendChild(igvCell);
                    
                    // Filter column with icons
                    const filterCell = document.createElement('td');
                    filterCell.className = 'col-filter';
                    const filterText = variant.Filter || 'N/A';
                    
                    if (filterText === 'PASS') {{
                        // Green tick for PASS
                        const icon = document.createElement('i');
                        icon.className = 'fas fa-check-circle';
                        icon.style.color = '#4CAF50';
                        icon.style.fontSize = '16px';
                        
                        // Center icon
                        const iconContainer = document.createElement('div');
                        iconContainer.className = 'center-icon';
                        iconContainer.appendChild(icon);
                        filterCell.appendChild(iconContainer);
                        
                        addTooltipToCell(filterCell, 'PASS');
                    }} else if (filterText === 'N/A' || filterText === '.') {{
                        // Create centered N/A container
                        const naContainer = document.createElement('div');
                        naContainer.className = 'center-na';
                        naContainer.textContent = 'N/A';
                        filterCell.appendChild(naContainer);
                    }} else {{
                        // Red X for any other filter value
                        const icon = document.createElement('i');
                        icon.className = 'fas fa-times-circle';
                        icon.style.color = '#E74C3C';
                        icon.style.fontSize = '16px';
                        
                        // Center icon
                        const iconContainer = document.createElement('div');
                        iconContainer.className = 'center-icon';
                        iconContainer.appendChild(icon);
                        filterCell.appendChild(iconContainer);
                        
                        addTooltipToCell(filterCell, filterText);
                    }}
                    
                    row.appendChild(filterCell);
                    
                    tbody.appendChild(row);
                }});
            }}
            
            // Render pagination controls
            function renderPagination() {{
                const pageCount = Math.ceil(filteredVariants.length / rowsPerPage);
                const pagination = document.getElementById('pagination');
                pagination.innerHTML = '';
                
                // Add Previous button
                const prevBtn = document.createElement('button');
                prevBtn.innerHTML = '<i class="fas fa-chevron-left"></i>';
                prevBtn.disabled = currentPage === 1;
                prevBtn.addEventListener('click', () => {{
                    if (currentPage > 1) {{
                        currentPage--;
                        renderTable(currentPage);
                        renderPagination();
                    }}
                }});
                pagination.appendChild(prevBtn);
                
                // Add page number buttons
                const maxButtons = 5;
                let startPage = Math.max(1, currentPage - Math.floor(maxButtons / 2));
                let endPage = Math.min(pageCount, startPage + maxButtons - 1);
                
                if (endPage - startPage + 1 < maxButtons && startPage > 1) {{
                    startPage = Math.max(1, endPage - maxButtons + 1);
                }}
                
                for (let i = startPage; i <= endPage; i++) {{
                    const pageButton = document.createElement('button');
                    pageButton.textContent = i;
                    pageButton.className = i === currentPage ? 'active' : '';
                    pageButton.addEventListener('click', () => {{
                        currentPage = i;
                        renderTable(currentPage);
                        renderPagination();
                    }});
                    pagination.appendChild(pageButton);
                }}
                
                // Add Next button
                const nextBtn = document.createElement('button');
                nextBtn.innerHTML = '<i class="fas fa-chevron-right"></i>';
                nextBtn.disabled = currentPage === pageCount;
                nextBtn.addEventListener('click', () => {{
                    if (currentPage < pageCount) {{
                        currentPage++;
                        renderTable(currentPage);
                        renderPagination();
                    }}
                }});
                pagination.appendChild(nextBtn);
            }}
            
            // Search functionality
            document.getElementById('searchInput').addEventListener('input', function(e) {{
                const searchTerm = e.target.value.toLowerCase();
                
                if (searchTerm === '') {{
                    filteredVariants = [...allVariants];
                }} else {{
                    filteredVariants = allVariants.filter(variant => {{
                        // Search in multiple fields
                        return (
                            (variant['Ref.Gene'] && variant['Ref.Gene'].toLowerCase().includes(searchTerm)) ||
                            (variant.CLNHGVS && variant.CLNHGVS.toLowerCase().includes(searchTerm)) ||
                            (variant.Chr && `${{variant.Chr}}:${{variant.Start}}`.toLowerCase().includes(searchTerm)) ||
                            (variant.avsnp151 && variant.avsnp151.toLowerCase().includes(searchTerm)) ||
                            (variant.Orpha && variant.Orpha.toLowerCase().includes(searchTerm)) ||
                            (variant.CLNDN && variant.CLNDN.toLowerCase().includes(searchTerm))
                        );
                    }});
                }}
                
                // Reset to first page and update display
                currentPage = 1;
                renderTable(currentPage);
                renderPagination();
            }});
            
            // Sorting functionality
            const tableHeaders = document.querySelectorAll('#variantsTable th');
            tableHeaders.forEach(header => {{
                header.addEventListener('click', function() {{
                    const sortBy = this.getAttribute('data-sort');
                    if (!sortBy) return;
                    
                    // Toggle sort direction
                    const currentDirection = this.getAttribute('data-direction') || 'asc';
                    const newDirection = currentDirection === 'asc' ? 'desc' : 'asc';
                    this.setAttribute('data-direction', newDirection);
                    
                    // Clear direction from other headers
                    tableHeaders.forEach(h => {{
                        if (h !== this) {{
                            h.removeAttribute('data-direction');
                        }}
                    }});
                    
                    // Sort the variants
                    filteredVariants.sort((a, b) => {{
                        let aValue, bValue;
                        
                        // Handle special sorting cases
                        switch(sortBy) {{
                            case 'rank':
                                aValue = a.rank;
                                bValue = b.rank;
                                break;
                                
                            case 'variant':
                                aValue = formatVariant(a.Chr, a.Start, a.Ref, a.Alt);
                                bValue = formatVariant(b.Chr, b.Start, b.Ref, b.Alt);
                                break;
                                
                            case 'gene':
                                aValue = a['Ref.Gene'] || '';
                                bValue = b['Ref.Gene'] || '';
                                break;
                                
                            case 'rs':
                                aValue = a.avsnp151 || '';
                                bValue = b.avsnp151 || '';
                                break;
                                
                            case 'effect':
                                aValue = getEffect(a) || '';
                                bValue = getEffect(b) || '';
                                break;
                                
                            case 'clnhgvs':
                                aValue = a.CLNHGVS || '';
                                bValue = b.CLNHGVS || '';
                                break;
                                
                            case 'hgvs':
                                // Sort by HGVS C as primary and HGVS P as secondary
                                aValue = (a['ANN[0].HGVS_C'] || '') + (a['ANN[0].HGVS_P'] || '');
                                bValue = (b['ANN[0].HGVS_C'] || '') + (b['ANN[0].HGVS_P'] || '');
                                break;
                                
                            case 'clinvar':
                                const aClinicvar = getClinvarValue(a);
                                const bClinicvar = getClinvarValue(b);
                                aValue = clinvarPriority[aClinicvar] || 999;
                                bValue = clinvarPriority[bClinicvar] || 999;
                                break;
                                
                            case 'acmg':
                                aValue = acmgPriority[a.ACMG] || 999;
                                bValue = acmgPriority[b.ACMG] || 999;
                                break;
                                
                            case 'acmgCodes':
                                aValue = a['ACMG Criteria'] || '';
                                bValue = b['ACMG Criteria'] || '';
                                break;
                                
                            case 'afGnomad':
                                // Handle scientific notation and numeric comparison
                                aValue = parseFloat(a['Freq_gnomAD_genome_ALL']) || 0;
                                bValue = parseFloat(b['Freq_gnomAD_genome_ALL']) || 0;
                                break;
                                
                            case 'af':
                                aValue = parseFloat(a.VAF) || 0;
                                bValue = parseFloat(b.VAF) || 0;
                                break;
                                
                            case 'phenotype':
                                aValue = formatPhenotype(a.Orpha) || '';
                                bValue = formatPhenotype(b.Orpha) || '';
                                break;
                                
                            case 'disease':
                                aValue = getDisease(a) || '';
                                bValue = getDisease(b) || '';
                                break;
                                
                            case 'database':
                                // Database column has the same value for all variants
                                // so we'll just sort by gene name as a fallback
                                aValue = a['Ref.Gene'] || '';
                                bValue = b['Ref.Gene'] || '';
                                break;
                                
                            case 'igv':
                                // Sort by IGV availability (top 300 variants first)
                                aValue = a.rank <= 300 ? 0 : 1;
                                bValue = b.rank <= 300 ? 0 : 1;
                                break;
                                
                            case 'filter':
                                aValue = a.Filter || '';
                                bValue = b.Filter || '';
                                break;
                                
                            default:
                                aValue = a[sortBy] || '';
                                bValue = b[sortBy] || '';
                        }}
                        
                        // Perform the comparison
                        if (typeof aValue === 'number' && typeof bValue === 'number') {{
                            return newDirection === 'asc' ? aValue - bValue : bValue - aValue;
                        }} else {{
                            const aStr = String(aValue).toLowerCase();
                            const bStr = String(bValue).toLowerCase();
                            return newDirection === 'asc' ? 
                                aStr.localeCompare(bStr) : 
                                bStr.localeCompare(aStr);
                        }}
                    }});
                    
                    // Re-render the table
                    renderTable(currentPage);
                }});
            }});
            
            // Initialize
            renderTable(currentPage);
            renderPagination();
        }});
    </script>
</body>
</html>"""
    
    # Write the HTML file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"Generated report: {output_path}")
    return True

def main():
    """Main function to run from the command line"""
    if len(sys.argv) < 2:
        print("Usage: python lctmbrep.py <variants_tsv_file> [tmb_report_file] [tumor_purity_file] [stats_file] [output_html_file]")
        return
    
    variants_file = sys.argv[1]
    
    # Check if TMB report file is provided
    tmb_report_file = None
    if len(sys.argv) > 2:
        tmb_report_file = sys.argv[2]
    
    # Check if tumor purity file is provided
    tumor_purity_file = None
    if len(sys.argv) > 3:
        tumor_purity_file = sys.argv[3]
    
    # Check if stats file is provided
    stats_file = None
    if len(sys.argv) > 4:
        stats_file = sys.argv[4]
    
    # Check if output file is provided
    output_file = None
    if len(sys.argv) > 5:
        output_file = sys.argv[5]
    
    if not os.path.exists(variants_file):
        print(f"Error: Variants file '{variants_file}' not found")
        return
    
    if tmb_report_file and not os.path.exists(tmb_report_file):
        print(f"Warning: TMB report file '{tmb_report_file}' not found")
    
    if tumor_purity_file and not os.path.exists(tumor_purity_file):
        print(f"Warning: Tumor purity file '{tumor_purity_file}' not found")
    
    if stats_file and not os.path.exists(stats_file):
        print(f"Warning: Stats file '{stats_file}' not found")
    
    result = generate_html(variants_file, tmb_report_file, tumor_purity_file, stats_file, output_file)
    if result:
        print("Report generation completed successfully!")
    else:
        print("Report generation failed.")

if __name__ == "__main__":
    main()
