import os
import re
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.sparse import issparse
import networkx as nx
import matplotlib.patches as patches
from matplotlib.patches import FancyArrowPatch
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from collections import defaultdict
import io
from PIL import Image

# -------------------------------  
# Define your positive/negative colors  
# -------------------------------  
default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']  
POS_COLOR = '#FF0056'   # for example
NEG_COLOR = '#008080'    # negative arrows  

def fmt_trunc_str(x, ndigits=2):
    s = str(x)
    if "." in s:
        whole, frac = s.split(".", 1)
        frac = (frac + "0"*ndigits)[:ndigits]
        return f"{whole}.{frac}"
    else:
        return f"{s}.{'0'*ndigits}"

# -------------------------------  
# Adjustable parameters  
# -------------------------------  
ICON_ZOOM                = 0.3
SHAPE_ICON_ZOOM          = 0.13
SHAPE_ARROW_LENGTH       = 1
GRID_SPACING             = 8
REGION_RADIUS            = 3
BASE_OFFSET              = 1.2
LABEL_PADDING            = 0.4
PARALLEL_PADDING         = 0.3
VERTICAL_LABEL_PAD       = 0.9
SCALE_INTRINSIC          = 2
INTRINSIC_ARROW_PADDING  = 0.8
CIRCLE_HEAD_RADIUS       = 0.15
CIRCLE_POINTER_SIZE      = CIRCLE_HEAD_RADIUS * 2
FIG_WIDTH_PER            = 5
FIG_HEIGHT               = 6

ION_CHANNEL_Y_OFFSETS = {1: -0.12}

SOURCE_LABEL_OFFSETS = {1:(0.0,0.0),2:(0.0,0.0),3:(0.0,0.0),4:(0.0,-0.5)}
TARGET_LABEL_OFFSETS = {1:(0.0,0.3),2:(0.0,0.0),3:(0.0,0.0),4:(0.5,-0.3)}
SELF_LABEL_OFFSETS   = {1:(0.0,0.0),2:(0.0,0.0),3:(0.0,0.0),4:(0.0,0.0)}

RAW_TITLE_MAP = {'covariate1':'baseline','covariate2':'sustained','covariate3':'transient'}

def draw_self_loop_with_circle(ax, position, loop_radius=0.4, color=NEG_COLOR, lw=2,
                               offset=(0,0), arrow_angle_deg=270, pointer_size=CIRCLE_POINTER_SIZE):
    cx, cy = position[0]+offset[0], position[1]+offset[1]
    loop = patches.Circle((cx,cy), loop_radius, edgecolor=color, facecolor='none', lw=lw, zorder=12)
    ax.add_patch(loop)
    angle = np.deg2rad(arrow_angle_deg)
    bx = cx+loop_radius*np.cos(angle)
    by = cy+loop_radius*np.sin(angle)
    head = patches.Circle((bx,by), pointer_size/2, edgecolor=color, facecolor=color, zorder=13)
    ax.add_patch(head)
    return loop, head

def process_file(mat_file_path):
    import re
    import numpy as np
    import scipy.io as sio
    from scipy.sparse import issparse

    # tolerant patterns (case-insensitive, allow spaces, optional 3rd index in T)
    pat_ex    = re.compile(r'(?i)\b(A|AN)\{(\d+)\}\(\s*(\d+)\s*,\s*(\d+)\s*\)')
    pat_shape = re.compile(r'(?i)^covariate\s+(\d+)\s*:\s*T\(\s*(\d+)\s*,\s*(\d+)(?:\s*,\s*\d+)?\s*\)$')

    def to_float(v):
        try:
            if issparse(v):
                return float(v.toarray()[0, 0])
        except Exception:
            pass
        try:
            return float(np.asarray(v).squeeze())
        except Exception:
            return float(v)

    def to_str(v):
        if isinstance(v, np.ndarray):
            return ''.join(v.astype(str)).strip()
        return str(v).strip()

    mat = sio.loadmat(mat_file_path, squeeze_me=True, struct_as_record=False)
    if 'Results' not in mat:
        raise KeyError(f"'Results' not found in {mat_file_path}")
    res = mat['Results']
    if not hasattr(res, 'PEB_thresholded'):
        raise KeyError(f"'Results.PEB_thresholded' not found in {mat_file_path}")
    peb = res.PEB_thresholded

    # flatten Ep/Pnames robustly for both single-struct and struct-array
    Ep_list, Pnames_list = [], []
    if isinstance(peb, np.ndarray):  # struct-array
        for elem in np.ravel(peb):
            Ep_list.append(to_float(getattr(elem, 'Ep')))
            Pnames_list.append(to_str(getattr(elem, 'Pnames')))
    else:  # single struct with array fields
        Ep_field = np.ravel(getattr(peb, 'Ep'))
        Pn_field = np.ravel(getattr(peb, 'Pnames'))
        Ep_list  = [to_float(x) for x in Ep_field]
        Pnames_list = [to_str(x) for x in Pn_field]

    Ep     = np.asarray(Ep_list, dtype=float)
    Pnames = list(Pnames_list)

    Np = len(Pnames)
    if Np == 0:
        raise ValueError("No parameter names found in PEB_thresholded.")
    nseg = len(Ep) // Np if len(Ep) % Np == 0 else 1

    segments = []
    for s in range(nseg):
        Ep_s = Ep[s*Np:(s+1)*Np] if nseg > 1 else Ep

        cov_ex    = {}
        cov_int   = {}
        cov_shape = {}

        for i, pn in enumerate(Pnames):
            val = Ep_s[i]

            # Always record T(...) (ion-channel) entries even if val == 0
            m_sh = pat_shape.search(pn)
            if m_sh:
                cov_idx, area, shape_idx = map(int, m_sh.groups())
                key = f"Covariate {cov_idx}"
                cov_shape.setdefault(key, []).append((area, shape_idx, val))
                # don't 'continue' — we may still parse extrinsic/intrinsic below for nonzero vals

            # skip zeros for extrinsic/intrinsic (keep your old behavior)
            if val == 0:
                continue

            # normalize "Covariate X: <conn>"
            if ':' in pn:
                cov, conn = re.split(r'\s*:\s*', pn, maxsplit=1)
                cov = cov.strip(); conn = conn.strip()
            else:
                cov, conn = pn, ''

            # extrinsic: A/AN{ins}(tgt,src)
            m_ex = pat_ex.search(conn)
            if m_ex:
                typ, ins, tgt, src = m_ex.groups()
                src_i, tgt_i = int(src), int(tgt)
                cov_ex.setdefault(cov, []).append((src_i, tgt_i, f"{typ}{ins}", val))

            # intrinsic: H(...), tied to a region r as third index in label
            if re.match(r'(?i)^covariate\s+\d+\b', cov) and conn.startswith('H'):
                cov_int.setdefault(cov, []).append({"label": conn, "Ep": val})

        segments.append({
            'cov_extrinsic':    cov_ex,
            'covariate_groups': cov_int,
            'cov_shape':        cov_shape,
            'segment':          s+1
        })

    return segments


region_positions = {1:(0,GRID_SPACING), 2:(GRID_SPACING,GRID_SPACING),
                    3:(0,0), 4:(GRID_SPACING,0)}
intrinsic_base_positions = {1:(-1,0.5),2:(0,2.0),3:(1,1.0),4:(0,0)}
cell_nodes = {1:'SS',2:'SP',3:'II',4:'DP'}
node_images = {1:r'D:\signal_data\fog\pics\Picture6.png',
               2:r'D:\signal_data\fog\pics\Picture5.png',
               3:r'D:\signal_data\fog\pics\Picture7.png',
               4:r'D:\signal_data\fog\pics\Picture8.png'}
shape_images = {1:r'D:\signal_data\fog\pics\Picture10.png',
                2:r'D:\signal_data\fog\pics\Picture20.png',
                3:r'D:\signal_data\fog\pics\Picture30.png'}

def plot_connectivity_for_covariate(ax, cov, data):
    ax.set_aspect('equal','box'); ax.axis('off')
    # region circles
    for pos in region_positions.values():
        ax.add_patch(patches.Circle(pos, REGION_RADIUS, facecolor='none', edgecolor='none', zorder=1))
    # ion channel shapes + arrows
    cov_shapes = data.get('cov_shape', {}).get(cov, [])
    shapes_by_area = defaultdict(list)
    for area, shape_idx, ep_val in cov_shapes:
        shapes_by_area[area].append((shape_idx, ep_val))
    base_angles = {1:135,2:45,3:225,4:315}
    spread_deg = 30
    for area, items in shapes_by_area.items():
        items_sorted = sorted(items, key=lambda x: x[0])
        angles = (np.linspace(base_angles[area]-spread_deg/2, base_angles[area]+spread_deg/2, len(items_sorted))
                  if len(items_sorted)>1 else [base_angles[area]])
        if area in (1,3):
            angles = angles[::-1]
        cx, cy = region_positions[area]
        for (shape_idx, ep_val), angle in zip(items_sorted, angles):
            dx = np.cos(np.deg2rad(angle))*(REGION_RADIUS+0.5)
            dy = np.sin(np.deg2rad(angle))*(REGION_RADIUS+0.5)
            x0, y0 = cx+dx, cy+dy + ION_CHANNEL_Y_OFFSETS.get(shape_idx, 0)

            img_arr = plt.imread(shape_images[shape_idx])
            img_box = OffsetImage(img_arr, zoom=SHAPE_ICON_ZOOM)
            if ep_val == 0:
                img_box.set_alpha(0.8)
            ab = AnnotationBbox(img_box, (x0, y0), frameon=False, zorder=5)
            ax.add_artist(ab)

            # positive or negative arrow
            if ep_val > 0:
                col = POS_COLOR
                ax.annotate('', xy=(x0, y0+SHAPE_ARROW_LENGTH), xytext=(x0, y0),
                            arrowprops=dict(arrowstyle='-|>', color=col, lw=1.5),
                            zorder=6)
                ax.text(x0, y0+SHAPE_ARROW_LENGTH, fmt_trunc_str(ep_val,2),
                        color=col, fontsize=8, va='bottom', ha='center',
                        bbox=dict(facecolor='white', alpha=0.7), zorder=6)
            elif ep_val < 0:
                col = NEG_COLOR
                ax.annotate('', xy=(x0, y0-SHAPE_ARROW_LENGTH), xytext=(x0, y0),
                            arrowprops=dict(arrowstyle='-|>', color=col, lw=1.5),
                            zorder=6)
                ax.text(x0, y0-SHAPE_ARROW_LENGTH, fmt_trunc_str(ep_val,2),
                        color=col, fontsize=8, va='top', ha='center',
                        bbox=dict(facecolor='white', alpha=0.7), zorder=6)

    # extrinsic edges
    cov_ex = data['cov_extrinsic'].get(cov, [])
    grouped_ex = {}
    for src, tgt, lab, ep in cov_ex:
        key = tuple(sorted((src, tgt)))
        grouped_ex.setdefault(key, []).append((src, tgt, lab, ep))
    for conns in grouped_ex.values():
        for idx, (src, tgt, lab, ep) in enumerate(conns):
            base_idx = idx - (len(conns)-1)/2
            sign = 1 if src < tgt else -1
            off_fac = base_idx * sign
            P_src = np.array(region_positions[src])
            P_tgt = np.array(region_positions[tgt])
            vec = P_tgt - P_src
            L = np.linalg.norm(vec)
            unit = vec / L if L else np.zeros(2)
            perp = np.array([-unit[1], unit[0]])
            start = P_src + perp*BASE_OFFSET*off_fac + unit*REGION_RADIUS
            end   = P_tgt + perp*BASE_OFFSET*off_fac - unit*REGION_RADIUS
            color = POS_COLOR if ep>0 else NEG_COLOR
            ax.annotate('', xy=tuple(end), xytext=tuple(start),
                        arrowprops=dict(arrowstyle='->', lw=3, color=color, mutation_scale=20))
            mid      = (start + end)/2
            off_perp = perp * LABEL_PADDING * off_fac
            pad      = VERTICAL_LABEL_PAD if abs(unit[0])<0.1 else PARALLEL_PADDING
            off_para = unit * pad * base_idx if abs(unit[0])<abs(unit[1]) else np.zeros(2)
            lbl_pos  = mid + off_perp + off_para
            ax.text(lbl_pos[0], lbl_pos[1],
                    f"{lab} ({fmt_trunc_str(ep, 2)})",
                    color='black', fontsize=10,
                    ha='center', va='center',
                    bbox=dict(facecolor='white', alpha=0.7), zorder=4)

    # intrinsic edges + self-loops
    cov_int = data['covariate_groups'].get(cov, [])
    for roi in region_positions:
        G = nx.MultiDiGraph()
        for n in cell_nodes: G.add_node(n)
        for ent in cov_int:
            if f",{roi})" in ent['label']:
                d,s,r = map(int, re.findall(r"\((\d+),(\d+),(\d+)\)", ent['label'])[0])
                if r==roi: G.add_edge(s,d,Ep=ent['Ep'])
        cent = np.mean(list(intrinsic_base_positions.values()), axis=0)
        pos = {n: tuple(np.array(region_positions[roi]) + SCALE_INTRINSIC*(np.array(p)-cent))
               for n,p in intrinsic_base_positions.items()}

        nonself,selfe = [],[]
        for u,v,d in G.edges(data=True):
            (selfe if u==v else nonself).append((u,v,d))

        # non-self intrinsic
        grouped = {}
        for u,v,d in nonself:
            grouped.setdefault(tuple(sorted((u,v))), []).append((u,v,d))
        for edges in grouped.values():
            for i,(u,v,d) in enumerate(edges):
                rad = 0.3 + 0.05*(i-(len(edges)-1)/2)
                col = POS_COLOR if d['Ep']>0 else NEG_COLOR
                p0,p1 = np.array(pos[u]), np.array(pos[v])
                vec2, L2 = p1-p0, np.linalg.norm(p1-p0)
                if L2 > 2*INTRINSIC_ARROW_PADDING:
                    uv2 = vec2/L2
                    p0 += uv2*INTRINSIC_ARROW_PADDING
                    p1 -= uv2*INTRINSIC_ARROW_PADDING
                arrow = FancyArrowPatch(p0, p1,
                                        connectionstyle=f"arc3,rad={rad}",
                                        arrowstyle='-|>' if u!=3 else '-',
                                        mutation_scale=20,
                                        color=col, lw=2, zorder=12)
                ax.add_patch(arrow)
                if u==3:
                    ax.add_patch(patches.Circle(p1, CIRCLE_HEAD_RADIUS, facecolor=col, zorder=13))
                path = arrow.get_path().transformed(arrow.get_patch_transform())
                mp = path.interpolated(steps=100).vertices[50]
                offset = (np.array(SELF_LABEL_OFFSETS[u]) if u==v
                          else np.array(SOURCE_LABEL_OFFSETS.get(u,(0,0))) +
                               np.array(TARGET_LABEL_OFFSETS.get(v,(0,0))))
                mp += offset
                ax.text(mp[0], mp[1], fmt_trunc_str(d['Ep'],2),
                        color=col, fontsize=9,
                        ha='center', va='center',
                        bbox=dict(facecolor='white', alpha=0.7), zorder=13)

        # self-loop edges
        offsets = {1:(0,-0.7),2:(0.7,0),3:(0,0.7),4:(-0.7,0)}
        for u,_,d in selfe:
            p = pos[u]
            col = POS_COLOR if d['Ep']>0 else NEG_COLOR
            loop, head = draw_self_loop_with_circle(ax, p,
                                                    loop_radius=0.4,
                                                    color=col,
                                                    offset=offsets[u])
            lbl = (np.array([p[0]+offsets[u][0]*1.5,
                             p[1]+offsets[u][1]*1.5]) +
                   np.array(SELF_LABEL_OFFSETS.get(u,(0,0))))
            ax.text(lbl[0], lbl[1], fmt_trunc_str(d['Ep'],2),
                    color=col, fontsize=9,
                    ha='center', va='center',
                    bbox=dict(facecolor='white', alpha=0.7),
                    zorder=14)

        # draw cell‐type icons
        for n,coord in pos.items():
            try:
                im = plt.imread(node_images[n])
                icon = OffsetImage(im, zoom=ICON_ZOOM)
                ab = AnnotationBbox(icon, coord, frameon=False)
                ax.add_artist(ab)
            except FileNotFoundError:
                ax.scatter(*coord, s=200, color='gray', zorder=10)

    xs = [p[0] for p in region_positions.values()]
    ys = [p[1] for p in region_positions.values()]
    pad = REGION_RADIUS + 1
    ax.set_xlim(min(xs)-pad, max(xs)+pad)
    ax.set_ylim(min(ys)-pad, max(ys)+pad)
    ax.margins(0.1)

def save_all_formats(fig, base_name, base_outdir, dpi=300):
    """
    Save `fig` into outdir/pre_post/{RGB,CMYK,Grayscale} as:
      - base_name.svg
      - base_name.tiff (300 dpi)
      - base_name.jpg (300 dpi)
    SVG is identical across subfolders (vector, sRGB). TIFF/JPG are converted.
    """
    prepost_dir = os.path.join(base_outdir, 'pre_post')
    rgb_dir     = os.path.join(prepost_dir, 'RGB')
    cmyk_dir    = os.path.join(prepost_dir, 'CMYK')
    gray_dir    = os.path.join(prepost_dir, 'Grayscale')
    for d in (rgb_dir, cmyk_dir, gray_dir):
        os.makedirs(d, exist_ok=True)

    # 1) Save SVG into each subfolder (vector; color space conceptually sRGB)
    for d in (rgb_dir, cmyk_dir, gray_dir):
        svg_path = os.path.join(d, f"{base_name}.svg")
        fig.savefig(svg_path, format='svg', bbox_inches='tight', transparent=True)
        print(f"Saved SVG: {svg_path}")

    # 2) Render once to a high-res PNG buffer with white background (for raster modes)
    buf = io.BytesIO()
    fig.savefig(
        buf, format='png', dpi=dpi, bbox_inches='tight',
        facecolor='white', edgecolor='white', transparent=False
    )
    buf.seek(0)
    rgb_im = Image.open(buf).convert('RGB')

    # 3) RGB outputs
    rgb_im.save(os.path.join(rgb_dir,  f"{base_name}.tiff"),
                format='TIFF', dpi=(dpi, dpi), compression='tiff_lzw')
    rgb_im.save(os.path.join(rgb_dir,  f"{base_name}.jpg"),
                format='JPEG', quality=95, dpi=(dpi, dpi), optimize=True)

    # 4) CMYK outputs
    cmyk_im = rgb_im.convert('CMYK')
    cmyk_im.save(os.path.join(cmyk_dir, f"{base_name}.tiff"),
                 format='TIFF', dpi=(dpi, dpi), compression='tiff_lzw')
    cmyk_im.save(os.path.join(cmyk_dir, f"{base_name}.jpg"),
                 format='JPEG', quality=95, dpi=(dpi, dpi), optimize=True)

    # 5) Grayscale outputs
    gray_im = rgb_im.convert('L')
    gray_im.save(os.path.join(gray_dir, f"{base_name}.tiff"),
                 format='TIFF', dpi=(dpi, dpi), compression='tiff_lzw')
    gray_im.save(os.path.join(gray_dir, f"{base_name}.jpg"),
                 format='JPEG', quality=95, dpi=(dpi, dpi), optimize=True)

# -------------------------------  
# MAIN entry  
# -------------------------------
if __name__=='__main__':
    file_paths = [
        r'D:\signal_data\cleaning_git\to_plot\plot_cTBS_average_pre_post.mat',
        r'D:\signal_data\cleaning_git\to_plot\plot_iTBS_average_pre_post.mat',
        r'D:\signal_data\cleaning_git\to_plot\plot_Sham_average_pre_post.mat'
    ]
    all_segments = [process_file(fp) for fp in file_paths]
    covs = sorted(set().union(
        *[set(seg['cov_extrinsic'].keys()) |
          set(seg['covariate_groups'].keys()) |
          set(seg['cov_shape'].keys())
          for segs in all_segments for seg in segs]
    ))
    outdir = r'D:\signal_data\cleaning_git\to_plot'
    os.makedirs(outdir, exist_ok=True)

    n_files = len(file_paths)
    n_segs  = len(all_segments[0])
    labels  = [os.path.splitext(os.path.basename(fp))[0] for fp in file_paths]

    for si in range(n_segs):
        for cov in covs:
            fig, axs = plt.subplots(1, n_files,
                                    figsize=(FIG_WIDTH_PER*n_files, FIG_HEIGHT),
                                    squeeze=False)
            axs = axs.flatten()
            fig.patch.set_facecolor('none')
            for ax in axs: ax.set_facecolor('none')

            for fi in range(n_files):
                plot_connectivity_for_covariate(axs[fi], cov, all_segments[fi][si])
                axs[fi].set_title(labels[fi], fontsize=10)

            key = cov.replace(' ','').lower()
            title = RAW_TITLE_MAP.get(key, cov)
            fig.suptitle(f"Segment {si+1} – {title}", fontsize=14)
            fig.subplots_adjust(left=0.05, right=0.98, top=0.9,
                                bottom=0.05, wspace=0)

            base_name = f"segment{si+1}_{key}"
            save_all_formats(fig, base_name, outdir, dpi=300)


