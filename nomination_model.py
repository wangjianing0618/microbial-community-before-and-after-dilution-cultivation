
import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm

# ==============================
# Simulation Parameters
# ==============================
np.random.seed(42)
n_trials = 100
resource_counts = [1, 5, 10, 20]  # Four resource environments
m_j = 0.1
resource_initial = 100
time = np.linspace(0, 100, 100)  # 100 time points

# Species functional groups and parameter ranges
n_per_group = 20
groups = ["HighG_HighB", "LowG_HighB", "HighG_LowB", "LowG_LowB"]
growth_ranges = {
    "HighG_HighB": (1.3, 1.6),
    "LowG_HighB":  (0.8, 1.0),
    "HighG_LowB":  (1.3, 1.6),
    "LowG_LowB":   (0.8, 1.0)
}
group_fullnames = {
    "HighG_HighB": "High growth + High metabolic flexibility",
    "LowG_HighB":  "Low growth + High metabolic flexibility",
    "HighG_LowB":  "High growth + Low metabolic flexibility",
    "LowG_LowB":   "Low growth + Low metabolic flexibility"
}
group_colors = {
    "HighG_HighB": "#4E79A7",
    "LowG_HighB":  "#F28E2B",
    "HighG_LowB":  "#E15759",
    "LowG_LowB":   "#76B7B2"
}

# Store the relative abundance of each group in each trial under all resource environments
results = []

# Total iterations = 4 resource environments × 100 trials
total_iterations = len(resource_counts) * n_trials

for n_resources in tqdm(resource_counts, desc="Resource environments"):
    HH_breadth_min = int(np.ceil(n_resources / 2))
    HH_breadth_max = n_resources

    for trial in tqdm(range(1, n_trials + 1), desc=f"Trials @ {n_resources} res", leave=False):
        # —— Construct 80 species: 20 per group ——
        species_list = []
        for grp in groups:
            g_min, g_max = growth_ranges[grp]
            # Metabolic breadth range: HighB: [ceil(n/2), n], LowB: [1, floor(n/2)]
            if "HighB" in grp:
                b_min, b_max = HH_breadth_min, HH_breadth_max
            else:
                b_min, b_max = 1, HH_breadth_min if n_resources > 1 else 1

            for _ in range(n_per_group):
                growth = np.random.uniform(g_min, g_max)
                breadth = np.random.randint(b_min, b_max + 1) if b_min <= b_max else 1
                resources_can_use = np.random.choice(n_resources, size=breadth, replace=False).tolist()
                species_list.append({
                    "Group": grp,
                    "Growth_Rate": growth,
                    "Resources": resources_can_use
                })

        n_species = len(species_list)  # = 80

        # —— Construct Uptake matrix U[i,j] ——
        U = np.zeros((n_species, n_resources))
        for i, sp in enumerate(species_list):
            for j in sp["Resources"]:
                U[i, j] = np.random.uniform(0.5, 1.5)
        growth_rates = np.array([sp["Growth_Rate"] for sp in species_list])

        # —— Initial abundance: each species N0 = 1 ——
        N0 = np.ones(n_species)
        R0 = np.ones(n_resources) * resource_initial
        y0 = np.concatenate([N0, R0])

        # —— Define multi-resource ODE system ——
        def crm_multi(y, t):
            N = y[:n_species]
            R = y[n_species:]
            resource_sum = U.dot(R)
            dN = N * (growth_rates * resource_sum - m_j)
            dR = np.zeros(n_resources)
            for j in range(n_resources):
                dR[j] = - R[j] * np.sum(growth_rates * N * U[:, j])
            return np.concatenate([dN, dR])

        # —— Numerical integration of ODE ——
        sol = odeint(crm_multi, y0, time)
        N_end = sol[-1, :n_species]
        N_end[N_end < 0] = 0

        # —— Sum absolute abundance by group & convert to relative abundance (%) ——
        group_totals = {grp: 0.0 for grp in groups}
        for i, sp in enumerate(species_list):
            group_totals[sp["Group"]] += N_end[i]
        total_end = sum(group_totals.values()) or 1.0

        for grp in groups:
            results.append({
                "Environment_Resources": n_resources,
                "Trial": trial,
                "Group": grp,
                "Final_Rel_Abundance": group_totals[grp] / total_end * 100
            })

# Convert to DataFrame
df_results = pd.DataFrame(results)

# ==============================
# Four-panel visualization: stacked area plot
# ==============================
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16
plt.rcParams['legend.fontsize'] = 16

fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharey=True)
axes = axes.flatten()

for idx, n_resources in enumerate(resource_counts):
    ax = axes[idx]
    df_env = df_results[df_results["Environment_Resources"] == n_resources]
    pivot = df_env.pivot(index="Trial", columns="Group", values="Final_Rel_Abundance").sort_index()
    x = pivot.index.values
    groups_order = ["HighG_HighB", "LowG_HighB", "HighG_LowB", "LowG_LowB"]
    y_stack = [pivot[g].values for g in groups_order]
    colors = [group_colors[g] for g in groups_order]
    labels = [group_fullnames[g] for g in groups_order]

    ax.stackplot(x, y_stack, labels=labels, colors=colors, alpha=0.8)
    ax.set_title(f"{n_resources} Resources Environment", fontname="Times New Roman")
    ax.set_xlabel("Trial", fontname="Times New Roman")
    if idx % 2 == 0:
        ax.set_ylabel("Relative Abundance (%)", fontname="Times New Roman")
    ax.set_xlim(1, n_trials)
    ax.set_ylim(0, 100)

# Unified legend on the right side
legend_elements = [plt.Line2D([0], [0], color=group_colors[g], lw=10, label=group_fullnames[g]) 
                   for g in groups_order]
fig.legend(handles=legend_elements,
           loc='center left',
           bbox_to_anchor=(1.02, 0.5),
           title="Functional Groups",
           title_fontsize=12,
           frameon=False)

plt.subplots_adjust(right=0.8, wspace=0.2, hspace=0.3)
plt.tight_layout()

# Save the figure to file, dpi=300, remove extra white border
plt.savefig("functional_groups_stacked.png", dpi=300, bbox_inches="tight")

# Optionally display the figure
plt.show()
