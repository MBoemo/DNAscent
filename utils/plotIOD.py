import sys
import numpy as np
import matplotlib.pyplot as plt


def parse_iod_output(filepath):
    headers = {}
    sections = {}
    current_section = None
    current_data = []

    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('#'):
                parts = line[1:].split(' ', 1)
                if len(parts) == 2:
                    headers[parts[0]] = parts[1]
            elif line.startswith('>'):
                if current_section is not None:
                    sections[current_section] = current_data
                current_section = line[1:].rstrip(':')
                current_data = []
            elif line and current_section is not None:
                current_data.append(line)

    if current_section is not None:
        sections[current_section] = current_data

    return headers, sections


def main():
    if len(sys.argv) != 2:
        print("Usage: python plotIOD.py <iod_output_file>")
        sys.exit(1)

    filepath = sys.argv[1]
    headers, sections = parse_iod_output(filepath)

    median_iod = float(headers.get('MedianIOD', '0').split()[0])
    ci = headers.get('95ConfidenceInterval', None)

    data_behind = np.array([float(x) for x in sections.get('DataBehindDistances', [])])
    sim_behind = np.array([float(x) for x in sections.get('SimBehindDistances', [])])
    data_ahead = np.array([float(x) for x in sections.get('DataAheadDistances', [])])
    sim_ahead = np.array([float(x) for x in sections.get('SimAheadDistances', [])])

    landscape_raw = sections.get('Landscape', [])
    landscape_fr, landscape_iod, landscape_w = [], [], []
    for row in landscape_raw:
        parts = row.split('\t')
        landscape_fr.append(float(parts[0]))
        landscape_iod.append(float(parts[1]))
        landscape_w.append(float(parts[2]))
    landscape_fr = np.array(landscape_fr)
    landscape_iod = np.array(landscape_iod)
    landscape_w = np.array(landscape_w)

    fig, axes = plt.subplots(2, 3, figsize=(20, 9))

    # Panel 1: Objective landscape
    ax = axes[0,0]
    ax.plot(landscape_fr, landscape_w, 'o-', color='#2166ac', markersize=3, linewidth=1)
    ax.set_xscale('log')
    ax.set_xlabel('Firing rate (fr)')
    ax.set_ylabel('Wasserstein distance')
    ax.set_title('Objective landscape')

    # Panel 2: Objective vs IOD
    ax = axes[0,1]
    ax.plot(landscape_iod, landscape_w, 'o-', color='#b2182b', markersize=3, linewidth=1)
    ax.set_xlabel('Median IOD (kb)')
    ax.set_ylabel('Wasserstein distance')
    ax.set_title('Objective vs IOD')
    ax.axvline(median_iod, color='grey', linestyle='--', linewidth=0.8, label=f'Best IOD = {median_iod:.0f} kb')
    if ci:
        lo, hi = [float(x) for x in ci.split() if x.replace('.','',1).replace('-','',1).isdigit()]
        ax.axvspan(lo, hi, alpha=0.15, color='grey', label=f'95% CI [{lo:.0f}, {hi:.0f}] kb')
    ax.legend(fontsize=8)

    # Panel 3: Histogram — behind-fork distances
    ax = axes[0,2]
    if len(data_behind) > 0 or len(sim_behind) > 0:
        all_behind = np.concatenate([d for d in [data_behind, sim_behind] if len(d) > 0])
        bins_b = np.linspace(0, np.percentile(all_behind, 99), 40)
        if len(data_behind) > 0:
            ax.hist(data_behind, bins=bins_b, density=True, alpha=0.5, color='#2166ac', label='Data')
        if len(sim_behind) > 0:
            ax.hist(sim_behind, bins=bins_b, density=True, alpha=0.5, color='#b2182b', label='Simulated')
    ax.set_xlabel('Behind-fork distance (kb)')
    ax.set_ylabel('Density')
    ax.set_title('Behind-fork histograms')
    ax.legend(fontsize=8)

    # Panel 4: CDF comparison — behind-fork distances
    ax = axes[1,0]
    if len(data_behind) > 0:
        s = np.sort(data_behind)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#2166ac', linewidth=1.5, label='Data')
    if len(sim_behind) > 0:
        s = np.sort(sim_behind)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#b2182b', linewidth=1.5, label='Simulated')
    ax.set_xlabel('Behind-fork distance (kb)')
    ax.set_ylabel('Cumulative probability')
    ax.set_title('Behind-fork CDFs')
    ax.legend(fontsize=8)

    # Panel 5: CDF comparison — ahead-of-fork distances
    ax = axes[1,1]
    if len(data_ahead) > 0:
        s = np.sort(data_ahead)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#2166ac', linewidth=1.5, label='Data')
    if len(sim_ahead) > 0:
        s = np.sort(sim_ahead)
        ax.step(s, np.arange(1, len(s)+1)/len(s), where='post', color='#b2182b', linewidth=1.5, label='Simulated')
    ax.set_xlabel('Ahead-of-fork distance (kb)')
    ax.set_ylabel('Cumulative probability')
    ax.set_title('Ahead-of-fork CDFs')
    ax.legend(fontsize=8)

    # Panel 6: Histogram — ahead-of-fork distances
    ax = axes[1,2]
    if len(data_ahead) > 0 or len(sim_ahead) > 0:
        all_ahead = np.concatenate([d for d in [data_ahead, sim_ahead] if len(d) > 0])
        bins_a = np.linspace(0, np.percentile(all_ahead, 99), 40)
        if len(data_ahead) > 0:
            ax.hist(data_ahead, bins=bins_a, density=True, alpha=0.5, color='#2166ac', label='Data')
        if len(sim_ahead) > 0:
            ax.hist(sim_ahead, bins=bins_a, density=True, alpha=0.5, color='#b2182b', label='Simulated')
    ax.set_xlabel('Ahead-of-fork distance (kb)')
    ax.set_ylabel('Density')
    ax.set_title('Ahead-of-fork histograms')
    ax.legend(fontsize=8)

    plt.tight_layout()
    outpath = filepath.rsplit('.', 1)[0] + '_IODplots.pdf'
    plt.savefig(outpath)
    print(f"Saved to {outpath}")
    plt.show()


if __name__ == '__main__':
    main()

