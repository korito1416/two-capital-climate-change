import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

output_dir = "./figure/baseline"
output_file_1 = os.path.join(output_dir, "Climate_Sensitivity.png")
output_file_2 = os.path.join(output_dir, "Temperature Anomaly.png")
output_file_3 = os.path.join(output_dir, "Intensity Function.png")


if not os.path.exists(output_dir):
    try:
        os.makedirs(output_dir)
    except Exception as e:
        print(f"Error creating directory {output_dir}: {e}")
        sys.exit(1)

#### Figure 1

def plot_hist(graph_title):
    # Read the data
    theta_list = pd.read_csv('./data/model144.csv', header=None).to_numpy()[:, 0] / 1000.

    # Create the histogram
    plt.figure(figsize=(8, 6))
    plt.hist(theta_list * 1000, bins=np.arange(0.9, 2.85, 0.15), density=True, 
             color="red", edgecolor='grey', linewidth=1, alpha=0.5)

    # Customize the plot
    plt.xlim(0.8, 3)
    plt.ylim(0, 1.5)
    plt.xlabel("Climate Sensitivity", fontsize=20)
    plt.ylabel("Density", fontsize=20)

    plt.savefig(graph_title)
    plt.close()


# Example usage:
plot_hist(output_file_1)



##### Figure 2

def plot_damage(graph_title):
    Y = np.linspace(0., 3., 1000)
    y_underline = 1.5
    y_overline = 2 

    gamma_1 = 0.00017675
    gamma_2 = 2 * 0.0022    
    gamma_3_list = np.linspace(0, 1. / 3, 20)

    y_limits = [1.5, 2.0]
    colors = ['red', 'blue']
    
    def y_data(y_limit):
        LHS_ylimit = np.zeros((1000,20))
        
        i=0
        for gamma_3 in gamma_3_list:
            LHS_ylimitlower = gamma_1 * Y + gamma_2/2 * Y**2 # y<y_limit
            LHS_ylimitupper = gamma_1 * Y + gamma_2*y_overline * \
                (Y-y_limit) + (gamma_2+gamma_3)/2 * \
                (Y-y_limit)**2 + gamma_2/2 * y_limit**2
        
            LHS_ylimit[:,i] =LHS_ylimitlower*(Y<y_limit) + LHS_ylimitupper*(Y>y_limit)
            i = i+1
        return LHS_ylimit

    plt.figure(figsize=(16, 10))

    for y_limit, color in zip(y_limits, colors):
        damages = y_data(y_limit)
        damage_upper = np.max(np.exp(-damages), axis=1)
        damage_lower = np.min(np.exp(-damages), axis=1)
        mean_damage = np.mean(np.exp(-damages), axis=1)

        plt.fill_between(Y, damage_lower, damage_upper, color=color, alpha=0.3)
        plt.plot(Y, mean_damage, color='black')
        plt.axvline(x=y_limit, color='black', linestyle="--")

    plt.ylim(0.65, 1)
    plt.xlim(0, 3)
    plt.xlabel("Temperature Anomaly", fontsize=20)
    plt.ylabel("Proportional reduction in economic output", fontsize=20)
    plt.savefig(graph_title)
    plt.grid(True)
    plt.show()

# Example usage:
plot_damage(output_file_2)


#### Figure 3


def plot_intensity(graph_title):
    def J(Y, y_underline=1.5):
        r1 = 1.5
        r2 = 2.5
        return r1 * (np.exp(r2 / 2 * (Y - y_underline)**2) - 1) * (Y >= y_underline)

    Y = np.linspace(0., 3., 1000)
    J_values = J(Y)

    plt.figure(figsize=(16, 10))

    plt.plot(Y, J_values, color="black", label=r"$\mathcal{J}(y)$",linewidth=6)

    plt.xlim(1, 2)
    plt.ylim(-0.01, 0.6)
    plt.xlabel("Temperature anomaly (ᵒC)", fontsize=20)
    plt.ylabel("Intensity", fontsize=20)

    plt.savefig(graph_title)


plot_intensity(output_file_3)






