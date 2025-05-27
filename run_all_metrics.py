import os
import subprocess
import sys
import pandas as pd
from scipy.stats import kruskal
import scikit_posthocs as sp

# Root directory containing all alignment data
ensemble_root = r"D:\DIPLOM\OrtDB"
# Directory where this wrapper script resides (and where we will save results and images)
script_dir = os.path.dirname(os.path.abspath(__file__))
results_dir = os.path.join(script_dir, 'results')
os.makedirs(results_dir, exist_ok=True)
# Path to the metrics script (assumed in the same directory as this wrapper)
metrics_script = os.path.join(script_dir, "alignment_metrics.py")

# Directory for saving images
imgs_dir = os.path.join(script_dir, 'imgs')
os.makedirs(imgs_dir, exist_ok=True)

# Directories to skip entirely
skip_dirs = {"sequences", "sum_sequences"}

print(f"Output CSV and heatmaps will be saved in: {script_dir}")
print(f"Heatmaps directory: {imgs_dir}\n")

#  НОВОЕ ДЛЯ OrtDB
output_files = []

# Walk through ensemble_root
for root, dirs, files in os.walk(ensemble_root):
    # Relative path under ensemble_root
    rel = os.path.relpath(root, ensemble_root)
    parts = rel.split(os.sep)
    # At the top level, skip special dirs
    if rel == ".":
        dirs[:] = [d for d in dirs if d not in skip_dirs]
        continue

    # Skip any directory under skip_dirs
    if parts[0] in skip_dirs:
        dirs[:] = []
        continue

    # Identify .fna files
    fna_files = [f for f in files if f.lower().endswith('.fna')]
    if not fna_files:
        continue

    # Determine sequence type
    dir_name = os.path.basename(root)
    if 'protein' in dir_name.lower() and any('Protein' in f for f in files):
        seq_type = 'aa'
    else:
        seq_type = 'nt'

    # Build output path in the script directory
    out_name = f"{dir_name}.csv"
    out_path = os.path.join(results_dir, out_name)

    # НОВОЕ ДЛЯ OrtDB
    output_files.append((dir_name, out_path))

    # Prepare command
    cmd = [sys.executable, metrics_script,
           '-d', root,
           '-t', seq_type,
           '-o', out_path,
           '--imgdir', imgs_dir]
    print(f"Running metrics on {root} (type={seq_type}) -> {out_path}")
    subprocess.run(cmd, check=True)

print("\nAll metrics computed and saved.")

# Собираем все результаты в одну длинную таблицу
long = []

# ВЕРНУТЬ ДЛЯ ENSEMBL

# for root, dirs, files in os.walk(results_dir):
#     # Тут results_dir — папка, куда вы сохраняли CSV
#     for fname in files:
#         if not fname.endswith('.csv'):
#             continue
#         sample = os.path.splitext(fname)[0]
#         df = pd.read_csv(os.path.join(results_dir, fname), sep=';', index_col=0).T
#         # df: строки = алгоритмы, столбцы = ['SP','TC','CS','Entropy']
#         for algo in df.index:
#             row = df.loc[algo].to_dict()
#             row.update({'Sample': sample, 'Algorithm': algo})
#             long.append(row)


for sample, csv_path in output_files:
    # csv_path уже указывает на существующий файл
    df = pd.read_csv(csv_path, sep=';', index_col=0).T
    # df: строки = алгоритмы, столбцы = ['SP','TC','CS','Entropy']
    for algo in df.index:
        row = df.loc[algo].to_dict()
        row.update({'Sample': sample, 'Algorithm': algo})
        long.append(row)


long_df = pd.DataFrame(long)
# Папка для статистики
stats_dir = os.path.join(script_dir, 'results')
os.makedirs(stats_dir, exist_ok=True)

metrics = ['SP', 'TC', 'CS', 'Entropy']
algos = sorted(long_df['Algorithm'].unique())
samples = sorted(long_df['Sample'].unique())

for metric in metrics:
    # Широкая таблица
    wide = long_df.pivot(index='Algorithm', columns='Sample', values=metric)
    wide.to_csv(
    os.path.join(stats_dir, f"{metric}_scores.csv"),
    sep=';',
    float_format='%.4f',
    index=True,
    index_label='Algorithm')

    print(f"\n=== {metric}-scores ===")
    print(wide.to_string(float_format='%.4f'))

    # Тест Краскала–Уоллиса
    groups = [long_df.loc[long_df['Algorithm'] == algo, metric].values for algo in algos]
    stat, p = kruskal(*groups)
    with open(os.path.join(stats_dir, f"kruskal_{metric}.txt"), 'w') as fh:
        fh.write(f"Kruskal–Wallis {metric}: H={stat:.4f}, p={p:.4g}\n")
    print(f"\nKruskal–Wallis on {metric}: H={stat:.4f}, p={p:.4g}")

    # Пост-хок Dunn с Bonferroni
    posthoc = sp.posthoc_dunn(long_df, val_col=metric, group_col='Algorithm', p_adjust='bonferroni')
    posthoc.to_csv(
    os.path.join(stats_dir, f"posthoc_{metric}.csv"),
    sep=';',
    float_format='%.4f')
    
    print(f"Post-hoc Dunn saved to posthoc_{metric}.csv")
