import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import os

def generate_plots():
    # 1. Load the data we generated in Step 3
    csv_file = "cyp2d6_analysis_results.csv"
    if not os.path.exists(csv_file):
        print(f"[ERROR] {csv_file} not found. Please run analysis.py first.")
        return
        
    df = pd.read_csv(csv_file)

    # ============================================================
    # PLOT 1: Seaborn Heatmap (Delta Instability)
    # ============================================================
    print("[INFO] Generating Seaborn Heatmap...")
    
    # Pivot the data so Alleles are rows and we plot the Delta_Instability
    # We set 'Allele' as index, and create a fake column for the heatmap
    heatmap_data = df.set_index("Allele")[["Delta_Instability"]]
    
    plt.figure(figsize=(8, 6))
    # cmap="RdBu_r" means Red is Destabilizing (positive), Blue is Stabilizing (negative)
    sns.heatmap(heatmap_data, annot=True, cmap="RdBu_r", center=0, 
                fmt=".2f", linewidths=.5, cbar_kws={'label': 'Change in Instability Index'})
    
    plt.title("Impact of CYP2D6 nsSNPs on Protein Stability\n(Red = Destabilizing, Blue = Stabilizing)", fontsize=14, pad=15)
    plt.ylabel("CYP2D6 Allele", fontsize=12)
    plt.xlabel("") # Hide X label since it's just one column
    plt.tight_layout()
    
    # Save the plot
    plt.savefig("stability_heatmap.png", dpi=300)
    print("[INFO] Saved stability_heatmap.png")
    plt.close()

    # ============================================================
    # PLOT 2: Plotly Interactive Scatter Plot
    # ============================================================
    print("[INFO] Generating Plotly Interactive Scatter Plot...")
    
    # We will plot Delta_MW on X, Delta_GRAVY on Y
    fig = px.scatter(
        df, 
        x="Delta_MW", 
        y="Delta_GRAVY",
        color="Effect_on_Stability",
        text="Allele",
        hover_data=["Mutation", "rsID", "Delta_Instability"],
        title="Physicochemical Changes in CYP2D6 Mutants",
        labels={
            "Delta_MW": "Change in Molecular Weight (Da)",
            "Delta_GRAVY": "Change in Hydrophobicity (GRAVY)",
            "Effect_on_Stability": "Stability Effect"
        },
        color_discrete_map={"Destabilizing": "red", "Stabilizing": "blue"}
    )
    
    # Improve layout
    fig.update_traces(textposition='top center', marker=dict(size=12, line=dict(width=2, color='DarkSlateGrey')))
    fig.update_layout(
        title_font_size=18,
        legend_title_font_size=12,
        template="plotly_white"
    )
    
    # Save as interactive HTML file
    fig.write_html("interactive_physicochemical_plot.html")
    print("[INFO] Saved interactive_physicochemical_plot.html")
    
    print("\n[SUCCESS] All plots generated successfully!")

if __name__ == "__main__":
    generate_plots()