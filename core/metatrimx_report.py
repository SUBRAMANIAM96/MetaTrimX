#!/usr/bin/env python3

import os
import re
import csv
import json
import datetime
import glob

# ==============================================================================
#                     METATRIMX REPORT GENERATOR v8.2
# ==============================================================================
# Author: SUBRAMANIAM VIJAYAKUMAR
# ==============================================================================
#
# SCRIPT FUNCTIONALITY:
# • Aggregates processing logs from all samples into a master CSV file.
# • Cross-references sequencing data with TensorFlow classification results.
# • Generates an interactive HTML dashboard for quality control and visualization.
#
# ==============================================================================

# --- CONFIGURATION ---
BASE_DIR = os.environ.get("OUTPUT_BASE_DIR", ".")
LOG_DIR = os.path.join(BASE_DIR, "Logs")
HTML_FILE = os.path.join(BASE_DIR, "MetaTrimX_Interactive_Report.html")
MASTER_LOG_FILE = os.path.join(BASE_DIR, "MetaTrimX_Master_Log.csv")
CLUST_LOG = os.path.join(LOG_DIR, "clustering_log.log")

# Path to the Neural Network classification output
AI_LOG = os.path.join(BASE_DIR, "Clustering_Results", "AI_Classification.csv")

# --- PARSING LOGIC ---
def parse_sample_log(log_path):
    """
    Parses an individual sample log file to extract read counts at each pipeline stage.
    """
    stats = {
        "Raw_Input": 0,
        "Demux_Output": 0,
        "Merged_Output": 0,
        "Final_Valid": 0
    }
    
    if not os.path.exists(log_path):
        return stats

    with open(log_path, 'r') as f:
        content = f.read()

    # 1. Raw Input
    raw_match = re.search(r"Total read pairs processed:\s+([\d,]+)", content)
    if raw_match:
        stats["Raw_Input"] = int(raw_match.group(1).replace(',', ''))

    # 2. Demux Output
    demux_match = re.search(r"Pairs written \(passing filters\):\s+([\d,]+)", content)
    if demux_match:
        stats["Demux_Output"] = int(demux_match.group(1).replace(',', ''))

    # 3. Merged Output
    merge_match = re.search(r"\s+([\d]+)\s+Merged", content)
    if merge_match:
        stats["Merged_Output"] = int(merge_match.group(1))

    # 4. Final Output
    final_match = re.search(r"(\d+) sequences kept", content)
    if final_match:
        stats["Final_Valid"] = int(final_match.group(1))

    return stats

def parse_clustering_stats():
    """Parses global clustering statistics from VSEARCH logs."""
    stats = {"uniques": 0, "otus": 0}
    if os.path.exists(CLUST_LOG):
        with open(CLUST_LOG, 'r') as f:
            content = f.read()
            u_match = re.search(r"(\d+)\s+unique sequences", content, re.IGNORECASE)
            if u_match: stats["uniques"] = int(u_match.group(1))
            
            o_match = re.search(r"Clusters:\s+(\d+)", content, re.IGNORECASE)
            if o_match: stats["otus"] = int(o_match.group(1))
            
            if stats["otus"] == 0:
                 asv_match = re.search(r"(\d+)\s+sequences generated", content, re.IGNORECASE)
                 if asv_match: stats["otus"] = int(asv_match.group(1))
    return stats

def parse_ml_stats():
    """
    Parses the Machine Learning classification file to count verified biological sequences.
    """
    stats = {"verified": 0, "total": 0}
    if os.path.exists(AI_LOG):
        try:
            with open(AI_LOG, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    stats["total"] += 1
                    # Count sequences classified as 'Real_Biology' by the model
                    if row.get("Verdict") == "Real_Biology":
                        stats["verified"] += 1
        except:
            pass
    return stats

def generate_master_log():
    """
    Scans the Logs directory, parses all sample logs, and writes the Master CSV.
    """
    print("[SYSTEM] Generating Master Log...")
    
    log_files = glob.glob(os.path.join(LOG_DIR, "*.log"))
    data_rows = []

    for log_file in log_files:
        filename = os.path.basename(log_file)
        if filename in ["clustering_log.log", "Pipeline_Global.log"]:
            continue
            
        sample_id = filename.replace(".log", "")
        stats = parse_sample_log(log_file)
        
        raw = stats["Raw_Input"] if stats["Raw_Input"] > 0 else 1
        demux_pct = round((stats["Demux_Output"] / raw) * 100, 2)
        demux_base = stats["Demux_Output"] if stats["Demux_Output"] > 0 else 1
        merge_pct = round((stats["Merged_Output"] / demux_base) * 100, 2)
        final_pct_raw = round((stats["Final_Valid"] / raw) * 100, 4)

        row = {
            "Sample_ID": sample_id,
            "Raw_Input": stats["Raw_Input"],
            "Demux_Output": stats["Demux_Output"],
            "Demux_Efficiency_%": demux_pct,
            "Merged_Output": stats["Merged_Output"],
            "Merge_Rate_%": merge_pct,
            "Final_Valid": stats["Final_Valid"],
            "Overall_Retention_%": final_pct_raw
        }
        data_rows.append(row)

    data_rows.sort(key=lambda x: x["Sample_ID"])

    if data_rows:
        headers = list(data_rows[0].keys())
        with open(MASTER_LOG_FILE, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers)
            writer.writeheader()
            writer.writerows(data_rows)
        print(f"[SUCCESS] Master Log written to: {MASTER_LOG_FILE}")
    else:
        print("[WARNING] No sample logs found to parse.")

    return data_rows

# --- HTML TEMPLATE ---
HTML_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>MetaTrimX Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600&display=swap" rel="stylesheet">
    <style>
        body {{ font-family: 'Inter', sans-serif; background: #f4f6f9; padding: 20px; color: #333; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 40px; border-radius: 8px; box-shadow: 0 4px 12px rgba(0,0,0,0.05); }}
        
        header {{ border-bottom: 2px solid #eee; padding-bottom: 20px; margin-bottom: 30px; display: flex; justify-content: space-between; align-items: center; }}
        h1 {{ margin: 0; color: #2c3e50; }}
        .meta {{ color: #7f8c8d; font-size: 14px; }}
        
        .stats-grid {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin-bottom: 40px; }}
        .card {{ background: #f8f9fa; padding: 20px; border-radius: 8px; border: 1px solid #eee; text-align: center; }}
        .card h3 {{ margin: 0; font-size: 32px; color: #2c3e50; }}
        .card span {{ color: #7f8c8d; font-size: 13px; text-transform: uppercase; letter-spacing: 1px; }}

        /* Badge for Machine Learning verification stats */
        .ml-badge {{ 
            display: inline-block;
            margin-top: 5px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
            color: white; 
            padding: 3px 8px; 
            border-radius: 12px; 
            font-size: 11px; 
            font-weight: bold;
        }}

        table {{ width: 100%; border-collapse: collapse; margin-top: 20px; font-size: 14px; }}
        th {{ background: #34495e; color: white; padding: 12px 15px; text-align: left; font-weight: 600; }}
        td {{ padding: 12px 15px; border-bottom: 1px solid #ddd; }}
        tr:nth-child(even) {{ background: #f9f9f9; }}
        
        .status-badge {{ padding: 4px 8px; border-radius: 4px; font-weight: 600; font-size: 12px; }}
        .good {{ background: #d5f5e3; color: #27ae60; }}
        .warn {{ background: #fcf3cf; color: #f39c12; }}
        .bad  {{ background: #fadbd8; color: #c0392b; }}
        
        .chart-section {{ margin-top: 50px; }}
        .chart-container {{ position: relative; height: 450px; width: 100%; }}
        .download-btn {{ background: #3498db; color: white; padding: 10px 20px; text-decoration: none; border-radius: 4px; font-weight: 600; display: inline-block; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div>
                <h1>MetaTrimX Analysis Report</h1>
                <div class="meta">Run Date: {date}</div>
            </div>
            <a href="{csv_link}" class="download-btn">Download Master Log (CSV)</a>
        </header>

        <div class="stats-grid">
            <div class="card">
                <h3>{total_samples}</h3>
                <span>Samples Processed</span>
            </div>
            <div class="card">
                <h3>{total_uniques}</h3>
                <span>Unique Sequences</span>
            </div>
            <div class="card">
                <h3>{total_otus}</h3>
                <span>Final OTUs</span>
                <br>
                <span class="ml-badge">✨ {ml_verified} Neural Verified</span>
            </div>
        </div>

        <h3>Sample Attrition Table</h3>
        <table>
            <thead>
                <tr>
                    <th>Sample ID</th>
                    <th>Raw Input</th>
                    <th>Primer Match (Demux)</th>
                    <th>Merged</th>
                    <th>Final Valid</th>
                    <th>Pipeline Efficiency</th>
                </tr>
            </thead>
            <tbody>
                {table_rows}
            </tbody>
        </table>

        <div class="chart-section">
            <h3>Data Loss Visualization (Attrition)</h3>
            <p style="color:#7f8c8d; font-size:14px;">This chart visualizes read retention across the pipeline steps. <br>
            <span style="color:#bdc3c7">■ Raw</span> -> <span style="color:#f39c12">■ Demux</span> -> <span style="color:#3498db">■ Merged</span> -> <span style="color:#27ae60">■ Final</span></p>
            <div class="chart-container">
                <canvas id="attritionChart"></canvas>
            </div>
        </div>
    </div>

    <script>
        const ctx = document.getElementById('attritionChart').getContext('2d');
        const data = {{
            labels: {js_labels},
            datasets: [
                {{
                    label: 'Final Valid',
                    data: {js_final},
                    backgroundColor: '#27ae60', // Green
                    order: 1
                }},
                {{
                    label: 'Merged',
                    data: {js_merged},
                    backgroundColor: '#3498db', // Blue
                    order: 2
                }},
                {{
                    label: 'Primer Match (Demux)',
                    data: {js_demux},
                    backgroundColor: '#f39c12', // Orange
                    order: 3
                }},
                {{
                    label: 'Raw Input',
                    data: {js_raw},
                    backgroundColor: '#bdc3c7', // Gray
                    order: 4
                }}
            ]
        }};

        new Chart(ctx, {{
            type: 'bar',
            data: data,
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                interaction: {{
                    mode: 'index',
                    intersect: false,
                }},
                scales: {{
                    x: {{ stacked: false }}, 
                    y: {{ 
                        type: 'logarithmic',
                        beginAtZero: true,
                        title: {{ display: true, text: 'Read Count (Log Scale)' }} 
                    }}
                }},
                plugins: {{
                    tooltip: {{
                        callbacks: {{
                            label: function(context) {{
                                return context.dataset.label + ': ' + context.raw.toLocaleString();
                            }}
                        }}
                    }}
                }}
            }}
        }});
    </script>
</body>
</html>
"""

def generate_html_report(data_rows):
    print("[SYSTEM] Generating HTML Dashboard...")
    
    rows_html = ""
    js_labels = []
    js_raw = []
    js_demux = []
    js_merged = []
    js_final = []

    c_stats = parse_clustering_stats()
    
    # Retrieve Machine Learning Stats
    ml_stats = parse_ml_stats()
    ml_display = ml_stats["verified"] if ml_stats["verified"] > 0 else "N/A"

    for row in data_rows:
        retention = row["Overall_Retention_%"]
        if retention > 10:
            status_class = "good"
            status_text = "OPTIMAL"
        elif retention > 1:
            status_class = "warn"
            status_text = "LOW YIELD"
        else:
            status_class = "bad"
            status_text = "FAILED"

        rows_html += f"""
        <tr>
            <td><strong>{row['Sample_ID']}</strong></td>
            <td>{row['Raw_Input']:,}</td>
            <td>
                {row['Demux_Output']:,}<br>
                <span style="font-size:11px; color:#7f8c8d;">{row['Demux_Efficiency_%']}%</span>
            </td>
            <td>
                {row['Merged_Output']:,}<br>
                <span style="font-size:11px; color:#7f8c8d;">{row['Merge_Rate_%']}%</span>
            </td>
            <td>
                <strong>{row['Final_Valid']:,}</strong>
            </td>
            <td>
                <span class="status-badge {status_class}">{retention}% Total</span>
            </td>
        </tr>
        """
        
        js_labels.append(row['Sample_ID'])
        js_raw.append(row['Raw_Input'])
        js_demux.append(row['Demux_Output'])
        js_merged.append(row['Merged_Output'])
        js_final.append(row['Final_Valid'])

    html_content = HTML_TEMPLATE.format(
        date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        csv_link="MetaTrimX_Master_Log.csv",
        total_samples=len(data_rows),
        total_uniques=f"{c_stats['uniques']:,}",
        total_otus=c_stats['otus'],
        # Inject ML verified count
        ml_verified=ml_display,
        table_rows=rows_html,
        js_labels=json.dumps(js_labels),
        js_raw=json.dumps(js_raw),
        js_demux=json.dumps(js_demux),
        js_merged=json.dumps(js_merged),
        js_final=json.dumps(js_final)
    )

    with open(HTML_FILE, 'w') as f:
        f.write(html_content)
    print(f"[SUCCESS] Report Generated: {HTML_FILE}")

# --- MAIN ---
if __name__ == "__main__":
    data = generate_master_log()
    if data:
        generate_html_report(data)