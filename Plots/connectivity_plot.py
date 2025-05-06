#!/usr/bin/env python3
"""
connectivity_plot.py

Generate brain‐region connectivity figures from SPM DCM PEB .mat files.

Usage:
    python src/connectivity_plot.py \
      --mat-files file1.mat file2.mat \
      [--outdir figures] \
      [--node-images SS.png SP.png II.png DP.png] \
      [--shape-images T1.png T2.png T3.png]
"""
import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.sparse import issparse
import networkx as nx
from matplotlib.patches import Circle, FancyArrowPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

# -------------------------------
# Project paths
# -------------------------------
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = BASE_DIR / 'data'
NODE_IMAGES_DEFAULT = {
    1: DATA_DIR / 'node_images' / 'SS.png',
    2: DATA_DIR / 'node_images' / 'SP.png',
    3: DATA_DIR / 'node_images' / 'II.png',
    4: DATA_DIR / 'node_images' / 'DP.png',
}
SHAPE_IMAGES_DEFAULT = {
    1: DATA_DIR / 'shape_images' / 'T1.png',
    2: DATA_DIR / 'shape_images' / 'T2.png',
    3: DATA_DIR / 'shape_images' / 'T3.png',
}

# -------------------------------
# Adjustable parameters
# -------------------------------
ICON_ZOOM               = 0.35
SHAPE_ICON_ZOOM         = 0.1
SHAPE_ARROW_LENGTH      = 1
GRID_SPACING            = 8
REGION_RADIUS           = 3
BASE_OFFSET             = 1.2
LABEL_PADDING           = 0.7
PARALLEL_PADDING        = 0.3
VERTICAL_LABEL_PAD      = 0.9
SCALE_INTRINSIC         = 2.1
INTRINSIC_ARROW_PADDING = 0.8
CIRCLE_HEAD_RADIUS      = 0.18
CIRCLE_POINTER_SIZE     = CIRCLE_HEAD_RADIUS * 2
FIG_WIDTH_PER           = 5
FIG_HEIGHT              = 6

ION_CHANNEL_Y_OFFSETS = { 1: -0.12 }
SELF_LOOP_OFFSETS = {
    1: (-0.4, 0.7),
    2: (0.7,  0.0),
    3: (0.0,  0.7),
    4: (0.8,  0.0),
}
SELF_LOOP_LABEL_OFFSETS = {
    1: (-0.4, 1.2),
    2: (1.0,  0.2),
    3: (0.0,  1.2),
    4: (1.0,  0.2),
}
RAW_TITLE_MAP = {
    'covariate1': 'baseline',
    'covariate2': 'sustained',
    'covariate3': 'transient',
}

REGION_POSITIONS = {
    1: (0, GRID_SPACING),
    2: (GRID_SPACING, GRID_SPACING),
    3: (0, 0),
    4: (GRID_SPACING, 0),
}
INTRINSIC_BASE_POSITIONS = {
    1: (-1, 0.5),
    2: (0, 2.0),
    3: (1, 1.0),
    4: (0, 0),
}
CELL_NODES = {1: "SS", 2: "SP", 3: "II", 4: "DP"}

# -------------------------------
# Utility functions
# -------------------------------
def fmt_trunc_str(x, ndigits=2):
    s = str(x)
    if "." in s:
        whole, frac = s.split(".", 1)
        frac = (frac + "0"*ndigits)[:ndigits]
        return f"{whole}.{frac}"
    return f"{s}.{'0'*ndigits}"

def draw_self_loop(ax, position, loop_radius=0.4, color='#008080',
                   lw=2, offset=(0, 0), angle_deg=270):
    cx, cy = position[0] + offset[0], position[1] + offset[1]
    loop = Circle((cx, cy), loop_radius, edgecolor=color, facecolor='none',
                  lw=lw, zorder=12)
    ax.add_patch(loop)
    theta = np.deg2rad(angle_deg)
    bx = cx + loop_radius * np.cos(theta)
    by = cy + loop_radius * np.sin(theta)
    head = Circle((bx, by), CIRCLE_POINTER_SIZE/2,
                  edgecolor=color, facecolor=color, zorder=13)
    ax.add_patch(head)

def process_file(mat_file):
    mat = sio.loadmat(mat_file, squeeze_me=True, struct_as_record=False)
    peb = mat['Results'].PEB_thresholded

    def to_float(v):
        return float(v.toarray()[0,0]) if issparse(v) else float(v)
    def to_str(v):
        return ''.join(v.astype(str)).strip() if isinstance(v, np.ndarray) else str(v)

    Ep     = np.array([to_float(x) for x in peb.Ep.flatten()])
    Pnames = [to_str(x)    for x in peb.Pnames.flatten()]
    Np     = len(Pnames)
    nseg   = len(Ep)//Np if len(Ep)%Np==0 else 1

    pat_ex    = re.compile(r'(A|AN)\{(\d+)\}\((\d+),(\d+)\)')
    pat_shape = re.compile(r'Covariate (\d+):T\((\d+),(\d+)\)')

    segments = []
    for s in range(nseg):
        Ep_s = Ep[s*Np:(s+1)*Np]
        cov_ex, cov_int, cov_shape = {}, {}, {}
        for val, pn in zip(Ep_s, Pnames):
            if val!=0 and ':' in pn:
                cov, conn = pn.split(':',1)
                m = pat_ex.search(conn)
                if m:
                    typ, ins, tgt, src = m.groups()
                    cov_ex.setdefault(cov, []).append((int(src),int(tgt),typ+ins,val))
            if val!=0 and pn.startswith("Covariate"):
                cov, conn = pn.split(':',1)
                if conn.strip().startswith("H"):
                    cov_int.setdefault(cov,[]).append({"label":conn.strip(),"Ep":val})
            m2 = pat_shape.match(pn)
            if m2:
                cov_idx, area, shape_idx = map(int, m2.groups())
                key = f"Covariate {cov_idx}"
                cov_shape.setdefault(key,[]).append((area,shape_idx,val))
        segments.append({
            'cov_extrinsic': cov_ex,
            'covariate_groups': cov_int,
            'cov_shape': cov_shape,
            'segment': s+1
        })
    return segments

def plot_covariate(ax, cov, data, node_images, shape_images):
    ax.set_aspect('equal','box')
    ax.axis('off')

    # draw regions
    for pos in REGION_POSITIONS.values():
        ax.add_patch(Circle(pos, REGION_RADIUS, facecolor='none',
                             edgecolor='black', lw=0, zorder=1))

    # ion-channel shapes & arrows
    cov_shapes = data.get('cov_shape', {}).get(cov, [])
    shapes_by_area = {}
    for area, shape_idx, ep_val in cov_shapes:
        shapes_by_area.setdefault(area, []).append((shape_idx, ep_val))
    base_angles = {1:135,2:45,3:225,4:315}
    spread_deg = 30
    for area, items in shapes_by_area.items():
        items_sorted = sorted(items, key=lambda x:x[0])
        n = len(items_sorted)
        angles = (np.linspace(base_angles[area]-spread_deg/2,
                              base_angles[area]+spread_deg/2,n)
                  if n>1 else [base_angles[area]])
        if area in (1,3): angles=angles[::-1]
        cx, cy = REGION_POSITIONS[area]
        for (shape_idx, ep_val), angle in zip(items_sorted, angles):
            dx = np.cos(np.deg2rad(angle))*(REGION_RADIUS+0.5)
            dy = np.sin(np.deg2rad(angle))*(REGION_RADIUS+0.5)
            x0, y0 = cx+dx, cy+dy
            y0 += ION_CHANNEL_Y_OFFSETS.get(shape_idx,0)

            # draw icon
            try:
                arr = plt.imread(shape_images[shape_idx])
                box = OffsetImage(arr, zoom=SHAPE_ICON_ZOOM)
                ab  = AnnotationBbox(box, (x0,y0), frameon=False, zorder=5)
                ax.add_artist(ab)
            except:
                ax.scatter(x0,y0,s=100,color='gray',zorder=5)

            # draw arrow + label
            if ep_val!=0:
                col = 'red' if ep_val>0 else 'blue'
                dy_arrow = SHAPE_ARROW_LENGTH if ep_val>0 else -SHAPE_ARROW_LENGTH
                ax.annotate('', xy=(x0,y0+dy_arrow), xytext=(x0,y0),
                            arrowprops=dict(arrowstyle='-|>',color=col,lw=1.5),
                            zorder=6)
                ax.text(x0, y0+dy_arrow, fmt_trunc_str(ep_val,2),
                        color=col, fontsize=8, ha='center',
                        va='bottom' if ep_val>0 else 'top',
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'),
                        zorder=6)

    # extrinsic connections
    cov_ex = data['cov_extrinsic'].get(cov, [])
    grouped = {}
    for src,tgt,lab,ep in cov_ex:
        grouped.setdefault(tuple(sorted((src,tgt))),[]).append((src,tgt,lab,ep))
    for conns in grouped.values():
        for idx,(src,tgt,lab,ep) in enumerate(conns):
            base_idx = idx - (len(conns)-1)/2
            sign = 1 if src<tgt else -1
            fac = base_idx*sign
            P0,P1 = np.array(REGION_POSITIONS[src]), np.array(REGION_POSITIONS[tgt])
            vec = P1-P0; L=np.linalg.norm(vec)
            unit = vec/L if L else np.zeros(2)
            perp = np.array([-unit[1],unit[0]])
            start = P0 + perp*BASE_OFFSET*fac + unit*REGION_RADIUS
            end   = P1 + perp*BASE_OFFSET*fac - unit*REGION_RADIUS
            col = '#FF0056' if ep>0 else '#008080'
            ax.annotate('', xy=tuple(end), xytext=tuple(start),
                        arrowprops=dict(arrowstyle='->',lw=3,color=col,mutation_scale=20))
            mid = (start+end)/2
            off_perp = perp*LABEL_PADDING*fac
            pad = VERTICAL_LABEL_PAD if abs(unit[0])<0.1 else PARALLEL_PADDING
            off_para = unit*pad*base_idx if abs(unit[0])<abs(unit[1]) else np.zeros(2)
            lbl = mid + off_perp + off_para
            ax.text(lbl[0],lbl[1],f"{lab} ({fmt_trunc_str(ep,2)})",
                    color='black',fontsize=10,ha='center',va='center',zorder=4)

    # intrinsic & self-loops
    cov_int = data['covariate_groups'].get(cov, [])
    for roi in REGION_POSITIONS:
        G = nx.MultiDiGraph()
        for n in CELL_NODES: G.add_node(n)
        for ent in cov_int:
            if f",{roi})" in ent['label']:
                d,s,r = map(int,re.findall(r"\((\d+),(\d+),(\d+)\)",ent['label'])[0])
                if r==roi: G.add_edge(s,d,Ep=ent['Ep'])
        cent = np.mean(list(INTRINSIC_BASE_POSITIONS.values()),axis=0)
        pos = {n: tuple(np.array(REGION_POSITIONS[roi])+
                       SCALE_INTRINSIC*(np.array(p)-cent))
               for n,p in INTRINSIC_BASE_POSITIONS.items()}
        nonself,selfe = [],[]
        for u,v,d in G.edges(data=True):
            (selfe if u==v else nonself).append((u,v,d))
        # nonself
        grp = {}
        for u,v,d in nonself:
            grp.setdefault(tuple(sorted((u,v))),[]).append((u,v,d))
        for edges in grp.values():
            for i,(u,v,d) in enumerate(edges):
                rad = 0.3 + 0.05*(i - (len(edges)-1)/2)
                col = '#FF0056' if d['Ep']>0 else '#008080'
                p0,p1 = np.array(pos[u]),np.array(pos[v])
                vec2,L2 = p1-p0, np.linalg.norm(p1-p0)
                if L2>2*INTRINSIC_ARROW_PADDING:
                    uv2 = vec2/L2; p0 += uv2*INTRINSIC_ARROW_PADDING; p1 -= uv2*INTRINSIC_ARROW_PADDING
                arrow = FancyArrowPatch(p0,p1,
                                       connectionstyle=f"arc3,rad={rad}",
                                       arrowstyle='-|>' if u!=3 else '-',
                                       mutation_scale=20,color=col,lw=2,zorder=12)
                ax.add_patch(arrow)
                if u==3:
                    ax.add_patch(Circle(p1,CIRCLE_HEAD_RADIUS,facecolor=col,zorder=13))
                path = arrow.get_path().transformed(arrow.get_patch_transform())
                mp = path.interpolated(steps=100).vertices[50]
                ax.text(mp[0],mp[1],fmt_trunc_str(d['Ep'],2),
                        color=col,fontsize=9,ha='center',va='center',
                        bbox=dict(facecolor='white',alpha=0.7,edgecolor='none'),zorder=13)
        # self-loops
        for u,_,d in selfe:
            p = pos[u]; off=SELF_LOOP_OFFSETS.get(u,(0,0))
            col = '#FF0056' if d['Ep']>0 else '#008080'
            draw_self_loop(ax,p,loop_radius=0.4,color=col,offset=off)
            loff=SELF_LOOP_LABEL_OFFSETS.get(u,(0,0))
            lbl=(p[0]+loff[0],p[1]+loff[1])
            ax.text(lbl[0],lbl[1],fmt_trunc_str(d['Ep'],2),
                    color=col,fontsize=9,ha='center',va='center',
                    bbox=dict(facecolor='white',alpha=0.7,edgecolor='none'),
                    zorder=14)
        # node icons
        for n,coord in pos.items():
            try:
                img = plt.imread(node_images[n])
                icon = OffsetImage(img,zoom=ICON_ZOOM)
                ab = AnnotationBbox(icon,coord,frameon=False)
                ax.add_artist(ab)
            except:
                ax.scatter(*coord,s=200,color='gray',zorder=10)

    # finalize
    xs=[p[0] for p in REGION_POSITIONS.values()]
    ys=[p[1] for p in REGION_POSITIONS.values()]
    pad = REGION_RADIUS+1
    extra = SHAPE_ARROW_LENGTH + max(abs(v) for v in ION_CHANNEL_Y_OFFSETS.values())
    ax.set_xlim(min(xs)-pad-extra, max(xs)+pad+extra)
    ax.set_ylim(min(ys)-pad-extra, max(ys)+pad+extra)
    ax.margins(0.05)

def main():
    parser = argparse.ArgumentParser(
        description='Generate connectivity plots from PEB .mat files.'
    )
    parser.add_argument(
        '--mat-files', '-m', nargs='+', required=True,
        help='Paths to one or more .mat files'
    )
    parser.add_argument(
        '--outdir', '-o', default='figures',
        help='Directory to save figures (default: ./figures)'
    )
    parser.add_argument(
        '--node-images', nargs=4, metavar=('SS','SP','II','DP'),
        help='Override paths to node icon images'
    )
    parser.add_argument(
        '--shape-images', nargs='+',
        help='Override paths to shape icon images'
    )
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    node_imgs = (
        dict(enumerate(args.node_images, start=1))
        if args.node_images else NODE_IMAGES_DEFAULT
    )
    shape_imgs = (
        dict(enumerate(args.shape_images, start=1))
        if args.shape_images else SHAPE_IMAGES_DEFAULT
    )

    all_segments = [process_file(str(fp)) for fp in args.mat_files]
    covs = sorted({
        cov
        for segs in all_segments
        for seg in segs
        for cov in (*seg['cov_extrinsic'],
                    *seg['covariate_groups'],
                    *seg['cov_shape'])
    })

    for f_idx, segs in enumerate(all_segments, start=1):
        for seg in segs:
            fig, axes = plt.subplots(
                1, len(covs),
                figsize=(FIG_WIDTH_PER * len(covs), FIG_HEIGHT),
                squeeze=False
            )
            axs = axes.flatten()
            for ax in axs:
                ax.set_facecolor('none')

            for ax, cov in zip(axs, covs):
                plot_covariate(ax, cov, seg, node_imgs, shape_imgs)
                key = cov.replace(' ', '').lower()
                ax.set_title(RAW_TITLE_MAP.get(key, cov), fontsize=10)

            fig.suptitle(f"File {f_idx} – Segment {seg['segment']}", fontsize=14)
            fig.subplots_adjust(left=0.05, right=0.98,
                                top=0.9, bottom=0.05, wspace=0)
            fname = f"file{f_idx}_segment{seg['segment']}.svg"
            outpath = outdir / fname
            fig.savefig(outpath, format='svg', bbox_inches='tight', transparent=True)
            print(f"Saved: {outpath}")
    print("Done.")

if __name__ == '__main__':
    main()
