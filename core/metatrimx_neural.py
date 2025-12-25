#!/usr/bin/env python3
import os
import numpy as np
import math
from collections import Counter
import sys

# Suppress TensorFlow/Scikit-Learn verbose logs
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

try:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import OneClassSVM
    from sklearn.feature_extraction.text import CountVectorizer
    from sklearn.decomposition import PCA
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

# ==============================================================================
#                 METATRIMX HYBRID CLASSIFICATION MODULE
# ==============================================================================
# Author: SUBRAMANIAM VIJAYAKUMAR
# ==============================================================================
#
# MODULE DESCRIPTION:
# • Implements a Dual-Path Supervised Learning architecture for sequence validation.
# • Path A (Random Forest): Analyzes statistical features (Entropy, GC Content).
# • Path B (One-Class SVM): Analyzes sequence motifs via K-mer embedding and PCA.
# • Ensemble Voting: Combines probability scores from both models for final classification.
#
# ==============================================================================

# --- MODULE 1: STATISTICAL FEATURE EXTRACTION ---
def calculate_entropy(seq):
    """Calculates Shannon entropy of the nucleotide sequence."""
    cnt = Counter(seq)
    length = len(seq)
    entropy = 0
    for count in cnt.values():
        p = count / length
        entropy -= p * math.log2(p)
    return entropy

def extract_features(seq):
    """
    Generates a feature vector: [Sequence Length, GC Content, Shannon Entropy]
    """
    seq = seq.upper()
    gc = (seq.count('G') + seq.count('C')) / len(seq) if seq else 0
    return [len(seq), gc, calculate_entropy(seq)]

# --- MODULE 2: SEQUENCE EMBEDDING (K-MER ANALYSIS) ---
def extract_kmers(seq, k=5):
    """Decomposes sequence into overlapping k-mers for vectorization."""
    return [seq[i:i+k] for i in range(len(seq)-k+1)]

# --- ENSEMBLE EXECUTION ENGINE ---
def run_hybrid_engine(sequences, ids):
    """
    Trains and executes the hybrid model on the provided sequence dataset.
    Returns a dictionary mapping Sequence IDs to Confidence Scores (0.0 - 1.0).
    """
    if not SKLEARN_AVAILABLE:
        print("  [ERROR] Scikit-learn dependency missing. Classification skipped.", flush=True)
        return {}

    print(f"  [MODEL] Initializing Dual-Path Analysis on {len(sequences)} sequences...", flush=True)

    # ----------------PATH A: RANDOM FOREST CLASSIFIER----------------
    print(f"  [PATH A] Training Random Forest on Statistical Features...", flush=True)
    X_features = np.array([extract_features(s) for s in sequences])
    
    # Define Training Set (Assumes top 20% by abundance are high-confidence 'Real' proxies)
    top_limit = max(20, int(len(sequences) * 0.2))
    
    # Generate Synthetic Noise for Binary Classification Training
    import random
    noise_seqs = [''.join(random.choices(['A','C','G','T'], k=200)) for _ in range(len(sequences))]
    X_noise = np.array([extract_features(s) for s in noise_seqs])
    
    # Combine Real Proxy + Synthetic Noise
    X_train = np.concatenate((X_features[:top_limit], X_noise))
    y_train = np.concatenate((np.ones(top_limit), np.zeros(len(X_noise))))
    
    rf_model = RandomForestClassifier(n_estimators=100, max_depth=10)
    rf_model.fit(X_train, y_train)
    rf_probs = rf_model.predict_proba(X_features)[:, 1]

    # ----------------PATH B: ONE-CLASS SVM (ANOMALY DETECTION)----------------
    print(f"  [PATH B] Training One-Class SVM on K-mer Embeddings...", flush=True)
    vectorizer = CountVectorizer(analyzer=extract_kmers)
    X_kmers = vectorizer.fit_transform(sequences).toarray()
    
    # Dimensionality Reduction (PCA) for SVM efficiency
    pca = PCA(n_components=10)
    X_pca = pca.fit_transform(X_kmers)
    
    svm_model = OneClassSVM(kernel="rbf", gamma="scale", nu=0.05)
    svm_model.fit(X_pca)
    svm_dists = svm_model.decision_function(X_pca)
    # Normalize SVM distance to probability-like score (Sigmoid)
    svm_probs = 1 / (1 + np.exp(-svm_dists))

    # ----------------ENSEMBLE CONSENSUS----------------
    print(f"  [ENSEMBLE] Aggregating probability scores...", flush=True)
    results = {}
    for i, pid in enumerate(ids):
        # Average probability from both models
        final_score = (rf_probs[i] + svm_probs[i]) / 2
        results[pid] = final_score

    return results

def run_ai_classification(fasta_file):
    """
    Main entry point called by the Core script.
    Parses FASTA, runs the engine, and returns scores.
    """
    if not os.path.exists(fasta_file): return {}
    seqs = []
    ids = []
    
    # Standard FASTA Parsing
    with open(fasta_file, 'r') as f:
        header = ""
        seq = ""
        for line in f:
            if line.startswith(">"):
                if header: 
                    ids.append(header.replace(">", "").strip())
                    seqs.append(seq)
                header = line.strip()
                seq = ""
            else: 
                seq += line.strip()
        if header: 
            ids.append(header.replace(">", "").strip())
            seqs.append(seq)
    
    if not seqs: return {}
    
    return run_hybrid_engine(seqs, ids)