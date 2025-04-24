import json
from datetime import datetime

def generate_standalone_html(json_file_path, output_html_path, unstable_file_path):
    """
    Generate a standalone HTML report for MSI analysis with aesthetics that match the WES report.
    
    Args:
        json_file_path: Path to the JSON file containing MSI analysis data
        output_html_path: Path to save the HTML report
        unstable_file_path: Path to the file containing unstable sites data
    """
    # Read the MSI analysis data
    with open(json_file_path, 'r') as f:
        data = json.load(f)
    
    # Get thresholds from technical details
    min_coverage = data['technical_details']['coverage_thresholds']['minimum']
    optimal_coverage = data['technical_details']['coverage_thresholds']['optimal']
    
    # Get the current date for the report
    current_date = datetime.now().strftime("%B %d, %Y")
    
    # Determine MSI status class for styling
    msi_status_class = "badge-conflicting" if data['msi_status']['classification'] == 'MSI-H' else "badge-benign"
    
    # Read unstable sites from file
    unstable_sites = read_unstable_sites(unstable_file_path)
    
    # HTML template with embedded CSS that matches WES report style
    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MSI Report - {data['sample_id']}</title>
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0-beta3/css/all.min.css" rel="stylesheet">
    <style>
        /* Color Palette - Matched with WES report */
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
            
            /* MSI specific colors */
            --mss-color: #4CAF50;
            --msi-high-color: #E74C3C;
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
        
        /* Report Header - Matched with WES report */
        .report-header {{
            background-color: var(--deep-blue);
            color: var(--white-bg);
            padding: 8.5px 12px;
            border-radius: 10px 10px 0 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            display: flex;
            justify-content: center;
            align-items: center;
            width: 23%;
            max-width: 231px;
            margin-left: 0;
            margin-right: auto;
            height: 45px;
        }}
        
        .logo img {{
            height: 29px;
            width: auto;
        }}
        
        /* Info Card - Matched with WES report */
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
        
        /* Badge Styling - Matched with WES report */
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
        
        .badge-benign, .badge-pass {{
            background: linear-gradient(90deg, var(--lime) 0%, #8ab22f 100%);
        }}
        
        .badge-conflicting, .badge-warning {{
            background: linear-gradient(90deg, var(--orange) 0%, #D97706 100%);
        }}
        
        .badge-fail {{
            background: linear-gradient(90deg, var(--red) 0%, #B91C1C 100%);
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
        
        /* Stats and Metrics Section - Matched with WES report */
        .stats-chart-container {{
            background-color: var(--white-bg);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }}
        
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(240px, 1fr));
            gap: 15px;
        }}
        
        .metric-card {{
            background-color: var(--soft-gray);
            padding: 15px;
            border-radius: 8px;
            display: flex;
            flex-direction: column;
            gap: 10px;
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }}
        
        .metric-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 3px 6px rgba(0,0,0,0.1);
        }}
        
        .metric-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .metric-title {{
            font-size: 14px;
            font-weight: 600;
            color: var(--deep-blue);
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        
        .metric-title i {{
            color: var(--lime);
        }}
        
        .metric-value {{
            font-size: 1.8rem;
            font-weight: bold;
            color: var(--deep-blue);
        }}
        
        .metric-details {{
            font-size: 0.9rem;
            color: var(--medium-gray);
        }}
        
        /* Progress Bars */
        .progress-container {{
            margin-top: 10px;
        }}
        
        .progress-label {{
            display: flex;
            justify-content: space-between;
            margin-bottom: 5px;
            font-size: 12px;
        }}
        
        .progress-label span {{
            color: var(--medium-gray);
        }}
        
        .progress-bar {{
            height: 5px;
            background-color: var(--soft-gray);
            border-radius: 3px;
            overflow: hidden;
        }}
        
        .progress-fill {{
            height: 100%;
            transition: width 0.3s ease;
            border-radius: 3px;
        }}
        
        .progress-fill.good {{
            background: linear-gradient(90deg, var(--lime) 0%, #8ab22f 100%);
        }}
        
        .progress-fill.warning {{
            background: linear-gradient(90deg, var(--orange) 0%, #D97706 100%);
        }}
        
        .progress-fill.poor {{
            background: linear-gradient(90deg, var(--red) 0%, #B91C1C 100%);
        }}
        
        /* Quality Section - Matched with WES report */
        .coverage-section {{
            background-color: var(--white-bg);
            padding: 10px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin-bottom: 15px;
        }}
        
        .coverage-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 8px;
            border-bottom: 1px solid var(--light-border);
            padding-bottom: 5px;
        }}
        
        .coverage-header h3 {{
            color: var(--deep-blue);
            font-size: 14px;
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        
        .coverage-metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 8px;
        }}
        
        .metrics-card {{
            background-color: var(--soft-gray);
            border-radius: 8px;
            padding: 8px;
            box-shadow: 0 1px 2px rgba(0,0,0,0.05);
            transition: transform 0.2s, box-shadow 0.2s;
            height: 55px;
        }}
        
        .metrics-card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 3px 6px rgba(0,0,0,0.1);
        }}
        
        .metrics-title {{
            font-size: 11px;
            color: var(--medium-gray);
            margin-bottom: 3px;
            display: flex;
            align-items: center;
            gap: 5px;
        }}
        
        .metrics-title i {{
            color: var(--lime);
        }}
        
        .metrics-value {{
            font-size: 15px;
            font-weight: 600;
            color: var(--deep-blue);
        }}
        
        /* Tables - Matched with WES report */
        .variants-container {{
            background-color: var(--white-bg);
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            overflow: auto;
            margin-bottom: 20px;
            position: relative;
        }}
        
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
            table-layout: fixed;
        }}
        
        .data-table th,
        .data-table td {{
            padding: 12px 15px;
            text-align: left;
            position: relative;
            text-overflow: ellipsis;
            white-space: nowrap;
        }}
        
        .data-table th {{
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
        
        .data-table th::after {{
            content: '';
            position: absolute;
            right: 10px;
            top: 50%;
            transform: translateY(-50%);
            opacity: 0.5;
            transition: opacity 0.3s ease;
        }}
        
        .data-table th[data-direction="asc"]::after {{
            content: '\\f0d8';
            font-family: 'Font Awesome 6 Free';
            font-weight: 900;
            opacity: 1;
        }}
        
        .data-table th[data-direction="desc"]::after {{
            content: '\\f0d7';
            font-family: 'Font Awesome 6 Free';
            font-weight: 900;
            opacity: 1;
        }}
        
        .data-table th:hover {{
            background: linear-gradient(180deg, #2C5282 0%, #1A365D 100%);
        }}
        
        .data-table tbody tr:nth-child(even) {{
            background-color: var(--soft-gray);
        }}
        
        .data-table tbody tr:hover {{
            background-color: rgba(44, 166, 164, 0.05);
        }}
        
        .data-table td {{
            border-bottom: 1px solid var(--light-border);
            overflow: hidden;
        }}
        
        /* MSI Thresholds Info Box - With WES report styling */
        .msi-threshold-info {{
            margin-top: 1rem;
            padding: 0.5rem;
            background-color: var(--light-teal);
            border-radius: 0.5rem;
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            align-items: center;
        }}
        
        .threshold-label {{
            display: inline-flex;
            align-items: center;
            gap: 8px;
            font-size: 0.875rem;
            color: var(--dark-text);
        }}
        
        .threshold-indicator {{
            width: 12px;
            height: 12px;
            border-radius: 50%;
        }}
        
        .threshold-mss .threshold-indicator {{
            background-color: var(--lime);
        }}
        
        .threshold-msi .threshold-indicator {{
            background-color: var(--conflicting-red);
        }}
        
        /* Footer - Matched with WES report */
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
        
        /* Responsive Design - Enhanced from WES report */
        @media (max-width: 1200px) {{
            .container {{
                padding: 10px;
            }}
            
            .metrics-grid {{
                grid-template-columns: repeat(2, 1fr);
            }}
            
            .coverage-metrics-grid {{
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
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Report Header - Matched with WES report style -->
        <div class="report-header">
            <div class="logo">
                <img src="logo.png" alt="Company Logo" style="height: 29px; width: auto;">
            </div>
        </div>
        
        <!-- Info Card - Styled to match WES report -->
        <div class="info-card">
            <h3>Sample Information</h3>
            <div class="info-details">
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-dna"></i></div>
                    <div class="info-text"><strong>MSI Results</strong> • Sample ID: {data['sample_id']}</div>
                </div>
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-calendar-alt"></i></div>
                    <div class="info-text"><strong>Analysis Date:</strong> {data['analysis_date']} • <strong>Report Date:</strong> {current_date}</div>
                </div>
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-clipboard-list"></i></div>
                    <div class="info-text">This report contains the microsatellite instability status and analysis metrics for the sample.</div>
                </div>
                <div class="info-item">
                    <div class="info-icon"><i class="fas fa-sort-amount-down"></i></div>
                    <div class="info-text">Sites are analyzed for coverage quality and instability patterns to determine MSI status.</div>
                </div>
            </div>
            
            <!-- MSI Classification Thresholds - Redesigned to match WES report styling -->
            <div class="msi-threshold-info">
                <div class="threshold-label threshold-mss">
                    <div class="threshold-indicator"></div>
                    <span><strong>MSS:</strong> &lt; 20% instability</span>
                </div>
                <div class="threshold-label threshold-msi">
                    <div class="threshold-indicator"></div>
                    <span><strong>MSI-H:</strong> ≥ 20% instability</span>
                </div>
            </div>
        </div>

        <!-- Stats Container - Styled to match WES report -->
        <div class="stats-chart-container">
            <div class="metrics-grid">
                <!-- MSI Status Metric Card -->
                <div class="metric-card">
                    <div class="metric-header">
                        <h3 class="metric-title"><i class="fas fa-chart-pie"></i> MSI Status</h3>
                        <span class="badge {msi_status_class}">{data['msi_status']['classification']}</span>
                    </div>
                    <div class="metric-value">{data['msi_status']['instability_percentage']:.1f}%</div>
                    <div class="metric-details">Instability percentage across analyzed sites</div>
                    <div class="progress-container">
                        <div class="progress-label">
                            <span>0%</span>
                            <span>20%</span>
                            <span>100%</span>
                        </div>
                        <div class="progress-bar">
                            <div class="progress-fill {'warning' if data['msi_status']['classification'] == 'MSI-H' else 'good'}"
                                 style="width: {min(data['msi_status']['instability_percentage'], 100)}%">
                            </div>
                        </div>
                    </div>
                </div>

                <!-- Coverage Metrics Card -->
                <div class="metric-card">
                    <div class="metric-header">
                        <h3 class="metric-title"><i class="fas fa-layer-group"></i> Coverage Metrics</h3>
                    </div>
                    <div class="metric-value">{data['quality_metrics']['all_sites']['mean_coverage']:.1f}<small>x</small></div>
                    <div class="metric-details">Mean coverage across all microsatellite sites</div>
                    <div class="progress-container">
                        <div class="progress-label">
                            <span>{min_coverage}x min</span>
                            <span>{optimal_coverage}x optimal</span>
                        </div>
                        <div class="progress-bar">
                            <div class="progress-fill {'good' if data['quality_metrics']['all_sites']['mean_coverage'] >= min_coverage else 'poor'}"
                                 style="width: {min(data['quality_metrics']['all_sites']['mean_coverage'] / optimal_coverage * 100, 100)}%">
                            </div>
                        </div>
                    </div>
                </div>

                <!-- Unstable Sites Card -->
                <div class="metric-card">
                    <div class="metric-header">
                        <h3 class="metric-title"><i class="fas fa-exclamation-triangle"></i> Unstable Sites</h3>
                    </div>
                    <div class="metric-value">{data['msi_status']['unstable_sites']} <small>of {data['msi_status']['total_sites']}</small></div>
                    <div class="metric-details">Number of microsatellite sites showing instability</div>
                    <div class="progress-container">
                        <div class="progress-label">
                            <span>0 sites</span>
                            <span>{data['msi_status']['total_sites']} sites</span>
                        </div>
                        <div class="progress-bar">
                            <div class="progress-fill {'warning' if data['msi_status']['classification'] == 'MSI-H' else 'good'}"
                                 style="width: {(data['msi_status']['unstable_sites'] / data['msi_status']['total_sites'] * 100)}%">
                            </div>
                        </div>
                    </div>
                </div>
                
                <!-- Coverage Quality Card -->
                <div class="metric-card">
                    <div class="metric-header">
                        <h3 class="metric-title"><i class="fas fa-check-circle"></i> Coverage Quality</h3>
                    </div>
                    <div class="metric-value">
                        {(data['quality_metrics']['all_sites']['sites_above_min_coverage'] / 
                         data['quality_metrics']['all_sites']['total_sites'] * 100):.1f}%
                    </div>
                    <div class="metric-details">Percentage of sites with coverage above {min_coverage}x</div>
                    <div class="progress-container">
                        <div class="progress-label">
                            <span>0%</span>
                            <span>100%</span>
                        </div>
                        <div class="progress-bar">
                            <div class="progress-fill {'good' if (data['quality_metrics']['all_sites']['sites_above_min_coverage'] / data['quality_metrics']['all_sites']['total_sites']) > 0.9 else 'warning'}"
                                 style="width: {(data['quality_metrics']['all_sites']['sites_above_min_coverage'] / data['quality_metrics']['all_sites']['total_sites'] * 100)}%">
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>



        <!-- Unstable Sites Section - With WES report styling -->
        <div class="coverage-section">
            <div class="coverage-header">
                <h3><i class="fas fa-dna"></i> Unstable Sites</h3>
            </div>
            
            <div class="variants-container">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>Location</th>
                            <th>Repeat Unit</th>
                            <th>Times</th>
                            <th>Instability</th>
                            <th>Coverage</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
                        {generate_unstable_sites_table(unstable_sites, min_coverage)}
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Footer - Matched with WES report -->
        <div class="footer">
            <p></p>
            <p><small></small></p>
        </div>
    </div>

    <script>
        // Add sorting functionality to table headers
        document.querySelectorAll('.data-table th').forEach(header => {{
            header.addEventListener('click', function() {{
                const table = this.closest('table');
                const index = Array.from(this.parentNode.children).indexOf(this);
                const rows = Array.from(table.querySelectorAll('tbody tr'));
                const direction = this.getAttribute('data-direction') === 'asc' ? 'desc' : 'asc';
                
                // Remove sort indicators from all headers
                table.querySelectorAll('th').forEach(th => th.removeAttribute('data-direction'));
                
                // Set sort direction on current header
                this.setAttribute('data-direction', direction);
                
                // Sort rows
                rows.sort((a, b) => {{
                    const aValue = a.children[index].textContent.trim();
                    const bValue = b.children[index].textContent.trim();
                    
                    // Handle numeric sorting for certain columns
                    if (index === 3 || index === 4) {{
                        const aNum = parseFloat(aValue);
                        const bNum = parseFloat(bValue);
                        return direction === 'asc' ? aNum - bNum : bNum - aNum;
                    }}
                    
                    // Default string comparison
                    return direction === 'asc' ? 
                        aValue.localeCompare(bValue) : 
                        bValue.localeCompare(aValue);
                }});
                
                // Reorder rows
                rows.forEach(row => table.querySelector('tbody').appendChild(row));
            }});
        }});
    </script>
</body>
</html>
'''

    # Write the complete HTML to file
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

def read_unstable_sites(file_path):
    """
    Read unstable sites from the provided file.
    
    Args:
        file_path: Path to the file containing unstable sites data
    
    Returns:
        List of unstable sites data
    """
    unstable_sites = []
    with open(file_path, 'r') as f:
        # Skip header line
        headers = f.readline().strip().split('\t')
        
        # Read all lines and convert to dictionary
        for line in f:
            values = line.strip().split('\t')
            site = {}
            for i, header in enumerate(headers):
                if i < len(values):
                    site[header] = values[i]
            
            # Add the site to the list
            unstable_sites.append(site)
    
    return unstable_sites

def generate_unstable_sites_table(sites, min_coverage):
    """
    Generate HTML for unstable sites table rows with enhanced styling.
    
    Args:
        sites: List of unstable sites data
        min_coverage: Minimum coverage threshold
    
    Returns:
        HTML string for table rows
    """
    rows = []
    for site in sites:
        # Get values from site data
        chromosome = site['chromosome']
        location = site['location']
        repeat_unit = site['repeat_unit_bases']
        repeat_times = site['repeat_times']
        instability = float(site['pro_p'])
        coverage = int(site['CovReads'])
        threshold = float(site['threshold'])
        
        # Determine status based on coverage
        quality_status = 'PASS' if coverage >= min_coverage else 'FAIL'
        
        # Determine status class based on quality status
        status_class = 'badge-benign' if quality_status == 'PASS' else 'badge-fail'
        
        # Determine instability class
        instability_class = 'badge-warning' if instability >= 0.2 else 'badge-benign'
        
        # Generate table row with enhanced styling
        rows.append(f'''
            <tr>
                <td><span style="font-family: monospace; font-weight: 600;">{chromosome}:{location}</span></td>
                <td>
                    <span style="display: inline-block; padding: 2px 6px; background-color: var(--light-teal); 
                          border-radius: 4px; font-family: monospace; font-weight: 600; color: var(--deep-blue);">
                        {repeat_unit}
                    </span>
                </td>
                <td>x{repeat_times}</td>
                <td>
                    <span class="badge {instability_class}">
                        {instability*100:.1f}%
                    </span>
                </td>
                <td>
                    <span style="font-weight: {600 if coverage >= min_coverage else 400};">
                        {coverage}x
                    </span>
                </td>
                <td><span class="badge {status_class}">{quality_status}</span></td>
            </tr>
        ''')
    return '\n'.join(rows)

if __name__ == '__main__':
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python lcmsirep.py <sample_id>")
        sys.exit(1)
    
    sample_id = sys.argv[1]
    json_file = f"{sample_id}_MSI_analysis.json"
    unstable_file = f"{sample_id}_MSI_unstable"
    output_html = f"{sample_id}_MSI_report.html"
    
    try:
        generate_standalone_html(json_file, output_html, unstable_file)
        print(f"Report generated successfully: {output_html}")
    except Exception as e:
        print(f"Error generating report: {str(e)}")
        sys.exit(1)
