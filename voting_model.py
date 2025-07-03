
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from collections import Counter

# ======================================
# 1. Simulate 100 trials, calculate relative abundance and ID of 80 species per trial
# ======================================

np.random.seed(42)

n_trials = 100
n_resources = 10
resource_initial = 100
m_j = 0.1
time = np.linspace(0, 100, 100)

growth_ranges = {
    "HighG_HighB": (1.3, 1.6),
    "LowG_HighB":  (0.8, 1.0),
    "HighG_LowB":  (1.3, 1.6),
    "LowG_LowB":   (0.8, 1.0)
}
other_groups = ["LowG_HighB", "HighG_LowB", "LowG_LowB"]

n_HH = 20
n_per_other_group = 20
n_others_total = 3 * n_per_other_group

# Calculate metabolic breadth range
HH_breadth_min = int(np.ceil(n_resources / 2))
HH_breadth_max = n_resources
if n_resources > 1:
    LB_breadth_min = 1
    LB_breadth_max = HH_breadth_min
else:
    LB_breadth_min = 1
    LB_breadth_max = 1

# 1.1 Generate 20 fixed HH species
fixed_HH_species = []
for i in range(n_HH):
    g = np.random.uniform(*growth_ranges["HighG_HighB"])
    b = np.random.randint(HH_breadth_min, HH_breadth_max + 1)
    resources_can_use = np.random.choice(n_resources, size=b, replace=False).tolist()
    fixed_HH_species.append({
        "ID": f"HH_{i+1}",
        "Group": "HighG_HighB",
        "Growth_Rate": g,
        "Resources": resources_can_use
    })

# 1.2 Define ODE system
def crm_multi(y, t, U, growth_rates):
    n_species = growth_rates.size
    N = y[:n_species]
    R = y[n_species:]
    resource_sum = U.dot(R)
    dN = N * (growth_rates * resource_sum - m_j)
    dR = np.zeros(n_resources)
    for j in range(n_resources):
        uptake_sum = np.sum(growth_rates * N * U[:, j])
        dR[j] = - R[j] * uptake_sum
    return np.concatenate([dN, dR])

# Store: species ID, relative abundance, and trial number for each trial
all_ids = []
all_rel_abundances = []
all_trial_indices = []

for trial in range(1, n_trials + 1):
    # Generate 60 other random species
    other_species = []
    for grp in other_groups:
        g_min, g_max = growth_ranges[grp]
        is_HighB = ("HighB" in grp)
        for j in range(n_per_other_group):
            g = np.random.uniform(g_min, g_max)
            if is_HighB:
                b_min, b_max = HH_breadth_min, HH_breadth_max
            else:
                b_min, b_max = LB_breadth_min, LB_breadth_max
            b = np.random.randint(b_min, b_max + 1) if b_min <= b_max else 1
            resources_can_use = np.random.choice(n_resources, size=b, replace=False).tolist()
            short_grp = grp.replace("HighG_HighB","HH")                           .replace("LowG_HighB","LH")                           .replace("HighG_LowB","HL")                           .replace("LowG_LowB","LL")
            other_species.append({
                "ID": f"{short_grp}_t{trial}_{j+1}",
                "Group": grp,
                "Growth_Rate": g,
                "Resources": resources_can_use
            })

    # Merge 80 species
    species_list = fixed_HH_species + other_species
    n_species = len(species_list)

    # Construct U and growth_rates
    U = np.zeros((n_species, n_resources))
    for i, sp in enumerate(species_list):
        for r in sp["Resources"]:
            U[i, r] = np.random.uniform(0.5, 1.5)
    growth_rates = np.array([sp["Growth_Rate"] for sp in species_list])

    # Initial abundance and integration
    N0 = np.zeros(n_species)
    hh_initial_each = 0.95 / n_HH
    others_initial_each = 0.05 / n_others_total
    N0[:n_HH] = hh_initial_each
    N0[n_HH:] = others_initial_each
    R0 = np.ones(n_resources) * resource_initial

    y0 = np.concatenate([N0, R0])
    sol = odeint(crm_multi, y0, time, args=(U, growth_rates))

    # Final abundance and relative abundance
    N_end = sol[-1, :n_species].copy()
    N_end[N_end < 0] = 0
    total_abund = np.sum(N_end)
    rel_abundances = (N_end / total_abund) * 100 if total_abund > 0 else np.zeros_like(N_end)

    # Collect species ID, relative abundance, and trial number
    for i, sp in enumerate(species_list):
        all_ids.append(sp["ID"])
        all_rel_abundances.append(rel_abundances[i])
        all_trial_indices.append(trial)

# ======================================
# 2. Stacked scatter plot (rotated): 
#    - X: Trial
#    - Y: Relative abundance (%)
#    - Color: HH_1…HH_20 use color, others in gray
#    - Font: Times New Roman
# ======================================
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 12

# Color map
cmap = plt.cm.get_cmap('tab20', n_HH)
hh_colors = {f"HH_{i+1}": cmap(i) for i in range(n_HH)}
other_color = '#CCCCCC'

# Prepare coordinates and colors
x_vals = all_trial_indices
y_vals = all_rel_abundances
point_colors = [hh_colors[id] if id in hh_colors else other_color for id in all_ids]

plt.figure(figsize=(12, 6))
plt.scatter(x_vals, y_vals, c=point_colors, s=30, edgecolor='black')

plt.xlabel("Simulation (Trial)", fontname="Times New Roman")
plt.ylabel("Relative Abundance (%)", fontname="Times New Roman")
plt.title("Stacked Scatter: All 80 Species per Trial (Rotated)", fontname="Times New Roman")
plt.xlim(0.5, n_trials + 0.5)
plt.ylim(0, max(all_rel_abundances) * 1.05)

# Legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor=hh_colors[f"HH_{i+1}"], edgecolor='black', label=f"Species {i+1}")
    for i in range(n_HH)
]
legend_elements.append(Patch(facecolor=other_color, edgecolor='black', label="Other"))
plt.legend(handles=legend_elements, bbox_to_anchor=(1.02, 0.5), loc='center left', frameon=False)

plt.tight_layout()
plt.savefig("stacked_scatter_rotated.png", dpi=300, bbox_inches="tight")
plt.show()

# ======================================
# 3. New analysis: count how often each HH species becomes Top1
# ======================================
top_ids = []
index = 0
for trial in range(1, n_trials + 1):
    rel_segment = all_rel_abundances[index:index + 80]
    id_segment = all_ids[index:index + 80]
    index += 80
    top_idx = np.argmax(rel_segment)
    top_ids.append(id_segment[top_idx])

counter = Counter(top_ids)
hh_counts = [counter.get(f"HH_{i+1}", 0) for i in range(n_HH)]

# Bar plot
plt.figure(figsize=(10, 5))
labels = [f"Species {i+1}" for i in range(n_HH)]
colors = [cmap(i) for i in range(n_HH)]

plt.bar(labels, hh_counts, color=colors, edgecolor='black')
plt.xlabel("HH Species", fontname="Times New Roman")
plt.ylabel("Count as Top1", fontname="Times New Roman")
plt.title("Number of Times Each HH Species Became Top1 (100 Trials)", fontname="Times New Roman")
plt.xticks(rotation=45, fontname="Times New Roman")
plt.yticks(fontname="Times New Roman")

plt.tight_layout()
plt.savefig("hh_top1_counts.png", dpi=300, bbox_inches="tight")
plt.show()
