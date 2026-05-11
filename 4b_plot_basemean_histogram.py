# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

os.chdir('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_MeCP2_SMARCB1')

# %%
# Read the CSV file
df = pd.read_csv('DEA_NSC.csv')
df.head()

# %%
print(f"Minimum baseMean: {df['baseMean'].min()}")
print(f"Median baseMean: {df['baseMean'].median()}")
print(f"Maximum baseMean: {df['baseMean'].max()}")

# %%
lower_bound = 100.0
upper_bound = df['baseMean'].quantile(0.98)
print(f"Lower bound: {lower_bound}")
print(f"Upper bound: {upper_bound}")

# %%
print(f"Total number of rows: {df.shape[0]}")
print(f"Number of rows above {lower_bound}: {df[df['baseMean'] > lower_bound].shape[0]}")

# %%
# Filter out outliers
df_filtered = df[(df['baseMean'] > lower_bound) & (df['baseMean'] < upper_bound)]
# %%
# Create the plot with improved style
plt.figure(figsize=(10, 6))
sns.histplot(data=df_filtered, x='baseMean', bins=50)

# Customize the plot
plt.title('Distribution of Base Mean Values (Excluding Outliers)', fontsize=12)
plt.xlabel('Base Mean', fontsize=10)
plt.ylabel('Count', fontsize=10)

# Add grid for better readability
plt.grid(True, alpha=0.3)

# Use scientific notation for x-axis if numbers are large
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig('results/basemean_histogram.png', dpi=300, bbox_inches='tight') 
# %%
