#!/usr/bin/env python3
"""Normalize a GROMACS liquid mass-density profile by a bulk reference."""
import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--profile', type=Path, required=True, help='GROMACS XVG: z in nm, densities in kg/m^3')
    parser.add_argument('--bulk-report', type=Path, required=True, help='JSON with density_mean_kg_m3')
    parser.add_argument('--liquid-column', type=int, default=1, help='Liquid data column, zero-based (default: 1)')
    parser.add_argument('--output-prefix', type=Path, required=True)
    parser.add_argument('--label', default='Liquid')
    parser.add_argument('--window', default='Final 1 ns of NVT')
    args = parser.parse_args()
    import math
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    bulk = float(json.loads(args.bulk_report.read_text())['density_mean_kg_m3'])
    if not math.isfinite(bulk) or bulk <= 0:
        parser.error('Bulk density must be finite and positive')
    rows = [[float(v) for v in line.split()] for line in args.profile.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith(('#', '@'))]
    if not rows or args.liquid_column < 1 or any(len(r) <= args.liquid_column for r in rows):
        parser.error('Missing liquid density column or profile data')
    z = [r[0] for r in rows]
    rho = [r[args.liquid_column] for r in rows]
    normalized = [value / bulk for value in rho]
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    with args.output_prefix.with_suffix('.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['z_relative_to_cellulose_center_nm', 'liquid_density_kg_m3', 'liquid_density_g_cm3', 'rho_over_rho_bulk'])
        writer.writerows((x, r, r/1000, q) for x, r, q in zip(z, rho, normalized))
    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(z, normalized, color='#1766a4', linewidth=1.7, label=args.label)
    ax.axhline(1, color='#333333', linestyle='--', linewidth=1.2, label='Bulk reference = 1')
    ax.set(xlabel='z relative to cellulose center (nm)', ylabel=r'$\rho_{\mathrm{liquid}}(z)\,/\,\rho_{\mathrm{bulk}}$',
           xlim=(min(z), max(z)), ylim=(0, max(1.2, max(normalized)*1.1)))
    ax.set_title(f'{args.label} density at the cellulose interface\n{args.window}', fontsize=12)
    ax.grid(axis='y', alpha=0.18)
    ax.legend(frameon=False)
    fig.text(0.5, -0.025, f'Simulated bulk reference: {bulk/1000:.6f} g/cm³ (preliminary; some density drift remains)',
             ha='center', fontsize=9, color='#444444')
    for suffix in ('.png', '.pdf'):
        fig.savefig(args.output_prefix.with_suffix(suffix), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'Normalized {len(rows)} bins using bulk density {bulk:.6f} kg/m³')
    print(f'Saved {args.output_prefix}.csv, .png, and .pdf')


if __name__ == '__main__':
    main()
