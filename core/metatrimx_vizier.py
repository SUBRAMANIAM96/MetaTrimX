#!/usr/bin/env python3

import os
import sys
import glob
import re
import json
import datetime
import statistics

# ==============================================================================
#                  METATRIMX REPORTER: DASHBOARD GENERATOR
# ==============================================================================
# Author: SUBRAMANIAM VIJAYAKUMAR
# ==============================================================================
#
# SCRIPT FUNCTIONALITY:
# • Parses log files from the cleaning (Step 1) and clustering (Step 2) pipelines.
# • Aggregates statistics on read retention, demultiplexing, and OTU generation.
# • Renders interactive HTML dashboards for quality control and result visualization.
#
# ==============================================================================

# --- CONFIGURATION (SMART PATH DETECTION) ---
# 1. Check if the environment variable is set (from Run Script)
env_dir = os.environ.get("OUTPUT_BASE_DIR")
hardcoded_dir = "MetaTrimX_Output"

# 2. DECISION LOGIC: Where are the actual logs?
# If the Env Var folder has logs, use it.
# If not, but 'MetaTrimX_Output' has logs, SWITCH to that (Fixes the mismatch).
if env_dir and os.path.exists(os.path.join(env_dir, "Logs")) and glob.glob(os.path.join(env_dir, "Logs", "*.log")):
    BASE_DIR = env_dir
elif os.path.exists(os.path.join(hardcoded_dir, "Logs")) and glob.glob(os.path.join(hardcoded_dir, "Logs", "*.log")):
    BASE_DIR = hardcoded_dir
    print(f"\n\033[93m[SYSTEM] Notice: Logs found in '{hardcoded_dir}'. Redirecting report there.\033[0m")
else:
    # Fallback: create the folder if it doesn't exist
    BASE_DIR = env_dir if env_dir else hardcoded_dir
    if not os.path.exists(BASE_DIR):
        try: os.makedirs(BASE_DIR, exist_ok=True)
        except: pass

LOG_DIR     = os.path.join(BASE_DIR, "Logs")
LOG_DIR_2   = os.path.join(BASE_DIR, "06_OTUs")  # Clustering Logs Directory

# Save these files INSIDE the Base Directory
OUTPUT_HTML = os.path.join(BASE_DIR, "MetaTrimX_Pipeline_Dashboard.html")
OUTPUT_HTML_2 = os.path.join(BASE_DIR, "MetaTrimX_Clustering_Dashboard.html")
OUTPUT_JSON = os.path.join(BASE_DIR, "MetaTrimX_Stats.json")

# --- ANSI TERMINAL COLORS ---
C_CYAN   = "\033[96m"
C_GREEN  = "\033[92m"
C_YELLOW = "\033[93m"
C_RED    = "\033[91m"
C_RESET  = "\033[0m"

# ==============================================================================
#                       SECTION 1: LOG PARSING MODULE
# ==============================================================================

class LogParser:
    """
    Regex Parser to extract read counts from bioinformatics tool logs.
    Handles variations in Cutadapt and Vsearch output formats.
    """
    def __init__(self, log_dir):
        self.log_dir = log_dir
        self.data = {}
        self.global_stats = {
            "Total_Raw": 0,
            "Total_Final": 0,
            "Samples_Count": 0,
            "Failed_Samples": []
        }

    def extract_number(self, pattern, text):
        """Finds a number using regex and converts to int."""
        match = re.search(pattern, text)
        if match:
            # Remove commas and convert to int
            return int(match.group(1).replace(",", ""))
        return None

    def parse_file(self, filepath):
        sample_id = os.path.basename(filepath).split("_")[0]
        
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        stats = {
            "Raw": 0,
            "Demux": 0,
            "Primers": 0,
            "Adapters": 0,
            "Merged": 0,
            "Final": 0
        }

        # --- 1. RAW READS (From Start of Log) ---
        raw = self.extract_number(r"Total read pairs processed:\s+([\d,]+)", content)
        if raw: stats["Raw"] = raw

        # --- 2. DEMULTIPLEXING ---
        all_cutadapt_passes = re.findall(r"Pairs written \(passing filters\):\s+([\d,]+)", content)
        
        if len(all_cutadapt_passes) >= 1:
            stats["Demux"] = int(all_cutadapt_passes[0].replace(",", ""))
        
        # --- 3. PRIMER REMOVAL ---
        if len(all_cutadapt_passes) >= 2:
            stats["Primers"] = int(all_cutadapt_passes[1].replace(",", ""))
        else:
            stats["Primers"] = stats["Demux"] # Fallback

        # --- 4. ADAPTER REMOVAL ---
        if len(all_cutadapt_passes) >= 3:
            stats["Adapters"] = int(all_cutadapt_passes[2].replace(",", ""))
        else:
            stats["Adapters"] = stats["Primers"] # Fallback

        # --- 5. MERGING (VSEARCH) ---
        merged = self.extract_number(r"\s+([\d,]+)\s+Merged", content)
        if merged: stats["Merged"] = merged

        # --- 6. FILTERING (VSEARCH MAXEE) ---
        final = self.extract_number(r"(\d+) sequences kept", content)
        if final: stats["Final"] = final

        # --- DATA INTEGRITY CHECK ---
        if stats["Raw"] == 0 and stats["Demux"] > 0:
            stats["Raw"] = stats["Demux"]

        return sample_id, stats

    def run(self):
        # Check if directory exists before globbing
        if not os.path.exists(self.log_dir):
            return None
            
        log_files = glob.glob(os.path.join(self.log_dir, "*.log"))
        
        if not log_files:
            return None

        print(f"{C_CYAN}[SYSTEM] Parsing Step 1 (Cleaning) logs from: {self.log_dir}{C_RESET}")

        for log_f in log_files:
            sid, s_data = self.parse_file(log_f)
            self.data[sid] = s_data
            
            # Update Globals
            self.global_stats["Total_Raw"] += s_data["Raw"]
            self.global_stats["Total_Final"] += s_data["Final"]
            self.global_stats["Samples_Count"] += 1
            
            if s_data["Final"] < 100:
                self.global_stats["Failed_Samples"].append(sid)

        return self.data

class ClusteringParser:
    """ 
    Parses VSEARCH logs for Clustering Stats (Step 2) 
    Extracts: Combined Reads, Unique Sequences, Chimeras, and Final OTUs
    """
    def __init__(self, log_dir):
        self.log_dir = log_dir
        self.data = {}

    def extract_number(self, pattern, text):
        match = re.search(pattern, text)
        return int(match.group(1).replace(",", "")) if match else 0

    def run(self):
        if not os.path.exists(self.log_dir): return None
        log_files = glob.glob(os.path.join(self.log_dir, "*.log"))
        if not log_files: return None

        print(f"{C_CYAN}[SYSTEM] Parsing Step 2 (Clustering) logs from: {self.log_dir}{C_RESET}")
        for log_f in log_files:
            sid = os.path.basename(log_f).replace(".log", "")
            with open(log_f, 'r', errors='ignore') as f: content = f.read()
            
            stats = {"Combined": 0, "Uniques": 0, "Chimeras": 0, "OTUs": 0}
            
            # 1. Combined Reads (Input to Derep)
            m_comb = re.search(r"Dereplicating file.*?in (\d+) seqs", content, re.DOTALL)
            if m_comb: stats["Combined"] = int(m_comb.group(1))
            
            # 2. Uniques (Output of Derep)
            stats["Uniques"] = self.extract_number(r"(\d+) uniques written", content)
            
            # 3. Chimeras (Output of Uchime)
            stats["Chimeras"] = self.extract_number(r"Found (\d+) .* chimeras", content)
            
            # 4. OTUs (Output of Cluster)
            stats["OTUs"] = self.extract_number(r"Clusters: (\d+)", content)
            
            self.data[sid] = stats
        return self.data

# ==============================================================================
#                       SECTION 2: HTML REPORT GENERATOR
# ==============================================================================

class HTMLGenerator:
    """
    Generates a standalone HTML Dashboard using embedded CSS.
    Visualizes key performance indicators and per-sample statistics.
    """
    def __init__(self, data, globals, filename):
        self.data = data
        self.globals = globals
        self.filename = filename

    def get_color(self, percentage):
        if percentage >= 70: return "#00f0ff" # Cyan (High)
        if percentage >= 50: return "#f0e68c" # Yellow (Med)
        return "#ff2a6d" # Neon Red (Low)

    def generate_css(self):
        return """
        <style>
            :root {
                --bg-dark: #050505;
                --bg-panel: #101010;
                --text-main: #e0e0e0;
                --neon-cyan: #00f0ff;
                --neon-green: #39ff14;
                --neon-red: #ff2a6d;
                --glass: rgba(255, 255, 255, 0.05);
                --border: 1px solid #333;
            }
            body { 
                background-color: var(--bg-dark); 
                color: var(--text-main); 
                font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif;
                margin: 0; padding: 20px;
                background-image: radial-gradient(circle at 50% 50%, #111 0%, #000 100%);
            }
            /* HEADER */
            .header {
                text-align: center;
                padding: 20px;
                border-bottom: 2px solid var(--neon-cyan);
                margin-bottom: 30px;
                text-transform: uppercase;
                letter-spacing: 3px;
                text-shadow: 0 0 10px var(--neon-cyan);
            }
            h1 { margin: 0; font-size: 2.5rem; color: var(--neon-cyan); }
            .subtitle { font-size: 0.9rem; color: #888; margin-top: 5px; }

            /* KPIS */
            .kpi-grid {
                display: grid;
                grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
                gap: 20px;
                margin-bottom: 40px;
            }
            .kpi-card {
                background: var(--bg-panel);
                border: 1px solid var(--neon-cyan);
                padding: 20px;
                border-radius: 8px;
                text-align: center;
                box-shadow: 0 0 15px rgba(0, 240, 255, 0.1);
            }
            .kpi-val { font-size: 2rem; font-weight: bold; color: #fff; }
            .kpi-label { font-size: 0.8rem; color: #aaa; text-transform: uppercase; }

            /* MAIN TABLE */
            .table-container {
                background: var(--bg-panel);
                border-radius: 10px;
                padding: 20px;
                overflow-x: auto;
                box-shadow: 0 10px 30px rgba(0,0,0,0.5);
                border: var(--border);
            }
            table { width: 100%; border-collapse: collapse; min-width: 800px; }
            th { 
                text-align: left; padding: 15px; 
                color: var(--neon-cyan); 
                border-bottom: 1px solid #444; 
                text-transform: uppercase; font-size: 0.85rem;
            }
            td { 
                padding: 12px 15px; 
                border-bottom: 1px solid #222;
                font-family: 'Consolas', monospace;
                color: #ccc;
            }
            tr:hover td { background-color: var(--glass); color: #fff; }
            
            /* PROGRESS BARS */
            .bar-bg {
                background: #222;
                height: 6px;
                width: 100%;
                border-radius: 3px;
                margin-top: 8px;
                position: relative;
            }
            .bar-fill {
                height: 100%;
                border-radius: 3px;
                box-shadow: 0 0 5px currentColor;
                transition: width 1s ease-in-out;
            }
            .stat-cell { min-width: 120px; }
            .trend-good { color: var(--neon-green); }
            .trend-bad { color: var(--neon-red); }

            /* FOOTER */
            .footer {
                text-align: center;
                margin-top: 50px;
                font-size: 0.8rem;
                color: #555;
            }
        </style>
        """

    def render(self):
        # Calculate Global Retention
        if self.globals['Total_Raw'] > 0:
            total_retention = (self.globals['Total_Final'] / self.globals['Total_Raw']) * 100
        else:
            total_retention = 0

        html = f"""<!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>MetaTrimX Analytics Dashboard</title>
            {self.generate_css()}
        </head>
        <body>
            <div class="header">
                <h1>🧬 MetaTrimX Analytics Core</h1>
                <div class="subtitle">Bioinformatics Pipeline Analysis Report • {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}</div>
            </div>

            <div class="kpi-grid">
                <div class="kpi-card">
                    <div class="kpi-val">{self.globals['Samples_Count']}</div>
                    <div class="kpi-label">Samples Processed</div>
                </div>
                <div class="kpi-card">
                    <div class="kpi-val" style="color:var(--neon-green)">{self.globals['Total_Final']:,}</div>
                    <div class="kpi-label">Reads Rescued</div>
                </div>
                <div class="kpi-card">
                    <div class="kpi-val" style="color:{self.get_color(total_retention)}">
                        {total_retention:.1f}%
                    </div>
                    <div class="kpi-label">Global Retention</div>
                </div>
                <div class="kpi-card">
                    <div class="kpi-val" style="color:var(--neon-red)">{len(self.globals['Failed_Samples'])}</div>
                    <div class="kpi-label">Failed / Low Yield</div>
                </div>
            </div>

            <div class="table-container">
                <h3 style="color:white; margin-bottom:20px; border-left: 4px solid var(--neon-cyan); padding-left:10px;">ATTRITION ANALYTICS</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>Raw Input</th>
                            <th>Demux</th>
                            <th>Primers</th>
                            <th>Merged</th>
                            <th>Final (Filtered)</th>
                            <th>Efficiency</th>
                        </tr>
                    </thead>
                    <tbody>
        """

        # Sort samples by ID
        for sid in sorted(self.data.keys()):
            s = self.data[sid]
            
            # Avoid division by zero
            if s['Raw'] == 0:
                pct = 0
            else:
                pct = (s['Final'] / s['Raw']) * 100
            
            color = self.get_color(pct)
            
            # Generate row
            html += f"""
                <tr>
                    <td style="font-weight:bold; color:var(--neon-cyan)">{sid}</td>
                    <td>{s['Raw']:,}</td>
                    <td>{s['Demux']:,}</td>
                    <td>{s['Primers']:,}</td>
                    <td>{s['Merged']:,}</td>
                    <td style="font-weight:bold; color:#fff">{s['Final']:,}</td>
                    <td class="stat-cell">
                        <div style="display:flex; justify-content:space-between; font-size:0.8rem;">
                            <span>{pct:.1f}%</span>
                        </div>
                        <div class="bar-bg">
                            <div class="bar-fill" style="width: {pct}%; background-color: {color}; color: {color}"></div>
                        </div>
                    </td>
                </tr>
            """

        html += """
                    </tbody>
                </table>
            </div>

            <div class="footer">
                GENERATED BY METATRIMX ANALYTICS • MODE 1 COMPILER <br>
                Running on Python Core v3.0
            </div>
        </body>
        </html>
        """
        
        # Write File
        with open(self.filename, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"{C_GREEN}[SUCCESS] Cleaning Dashboard: {self.filename}{C_RESET}")

    def render_clustering(self):
        """ Generates Step 2 Clustering Report using the same CSS architecture """
        
        # Calculate Globals for KPI
        total_otus = sum(d['OTUs'] for d in self.data.values())
        total_chimeras = sum(d['Chimeras'] for d in self.data.values())
        
        html = f"""<!DOCTYPE html>
        <html lang="en">
        <head>
            <meta charset="UTF-8">
            <title>MetaTrimX Clustering Dashboard</title>
            {self.generate_css()}
        </head>
        <body>
            <div class="header">
                <h1>🦠 MetaTrimX Analytics Core</h1>
                <div class="subtitle">Clustering Analytics Report • {datetime.datetime.now().strftime("%Y-%m-%d %H:%M")}</div>
            </div>

            <div class="kpi-grid">
                <div class="kpi-card">
                    <div class="kpi-val">{len(self.data)}</div>
                    <div class="kpi-label">Samples Clustered</div>
                </div>
                <div class="kpi-card">
                    <div class="kpi-val" style="color:var(--neon-cyan)">{total_otus:,}</div>
                    <div class="kpi-label">Total OTUs Found</div>
                </div>
                <div class="kpi-card">
                    <div class="kpi-val" style="color:var(--neon-red)">{total_chimeras:,}</div>
                    <div class="kpi-label">Chimeras Removed</div>
                </div>
            </div>

            <div class="table-container">
                <h3 style="color:white; margin-bottom:20px; border-left: 4px solid var(--neon-cyan); padding-left:10px;">CLUSTERING STATS</h3>
                <table>
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>Combined Reads</th>
                            <th>Unique Seqs</th>
                            <th>Chimeras Detected</th>
                            <th>Final OTUs</th>
                            <th>Diversity %</th>
                        </tr>
                    </thead>
                    <tbody>
        """

        for sid in sorted(self.data.keys()):
            s = self.data[sid]
            # Diversity Score (OTUs / Uniques)
            if s['Uniques'] > 0:
                div = (s['OTUs'] / s['Uniques']) * 100
            else:
                div = 0
            
            html += f"""
                <tr>
                    <td style="font-weight:bold; color:var(--neon-cyan)">{sid}</td>
                    <td>{s['Combined']:,}</td>
                    <td>{s['Uniques']:,}</td>
                    <td style="color:var(--neon-red)">{s['Chimeras']:,}</td>
                    <td style="font-weight:bold; color:var(--neon-green); font-size:1.1em">{s['OTUs']:,}</td>
                    <td class="stat-cell">
                        <div class="bar-bg">
                            <div class="bar-fill" style="width: {div}%; background-color: var(--neon-cyan);"></div>
                        </div>
                    </td>
                </tr>
            """

        html += """
                    </tbody>
                </table>
            </div>

            <div class="footer">
                GENERATED BY METATRIMX ANALYTICS • MODE 2 CLUSTERING
            </div>
        </body>
        </html>
        """
        
        with open(self.filename, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"{C_GREEN}[SUCCESS] Clustering Dashboard: {self.filename}{C_RESET}")

# ==============================================================================
#                       SECTION 3: MAIN EXECUTION
# ==============================================================================

def main():
    print(f"{C_CYAN}")
    print("=========================================================")
    print("      📊 METATRIMX REPORTER: INITIALIZING... 📊")
    print("=========================================================")
    print(f"{C_RESET}")

    # --- REPORT 1: CLEANING (Step 1) ---
    if os.path.exists(LOG_DIR):
        parser = LogParser(LOG_DIR)
        data = parser.run()
        if data:
            # Generate the first report
            generator = HTMLGenerator(data, parser.global_stats, OUTPUT_HTML)
            generator.render()
            
            # Save JSON stats
            with open(OUTPUT_JSON, 'w') as f:
                json.dump(parser.data, f, indent=4)
            print(f"{C_GREEN}[SUCCESS] Stats exported to {OUTPUT_JSON}{C_RESET}")

    # --- REPORT 2: CLUSTERING (Step 2) [NEW] ---
    if os.path.exists(LOG_DIR_2):
        parser2 = ClusteringParser(LOG_DIR_2)
        data2 = parser2.run()
        if data2:
            # Generate the second report
            generator2 = HTMLGenerator(data2, {}, OUTPUT_HTML_2)
            generator2.render_clustering()

    print(f"{C_CYAN}")
    print("---------------------------------------------------------")
    print("      REPORTS COMPLETE. OPEN HTML FILES TO VIEW.")
    print("---------------------------------------------------------")
    print(f"{C_RESET}")

if __name__ == "__main__":
    main()