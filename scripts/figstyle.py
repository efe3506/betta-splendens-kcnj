"""Common Matplotlib style. Import this module before pyplot in every figure script."""
import os

import matplotlib

matplotlib.use('Agg')

matplotlib.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Liberation Sans', 'Nimbus Sans', 'DejaVu Sans'],
    'mathtext.fontset': 'custom',
    'mathtext.rm': 'Arial',
    'mathtext.it': 'Arial:italic',
    'mathtext.bf': 'Arial:bold',
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'savefig.dpi': 300,
    'savefig.facecolor': 'white',
    'savefig.bbox': 'tight',
    'figure.facecolor': 'white',
    'axes.linewidth': 0.9,
})

SHOW_TITLES = os.environ.get('FIG_TITLES', '0') == '1'


def title(ax, *args, **kwargs):
    if SHOW_TITLES:
        ax.set_title(*args, **kwargs)
