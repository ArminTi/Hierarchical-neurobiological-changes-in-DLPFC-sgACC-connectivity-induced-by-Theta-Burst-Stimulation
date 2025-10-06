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

# --- Added for extra saving options ---
import io
from PIL import Image
# --------------------------------------

# -----------------------------------------
# Arrow colors
# -----------------------------------------
POS_COLOR = '#FF0056'   # positive arrows
NEG_COLOR = '#008080'   # negative arrows

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
ICON_ZOOM               = 0.35    # size of the node icons
SHAPE_ICON_ZOOM         = 0.1     # zoom for ion-channel icons
SHAPE_ARROW_LENGTH      = 0.8     # length of ion-channel arrows
GRID_SPACING            = 10      # distance between region centers
REGION_RADIUS           = 3       # radius of each brain-area circle
BASE_OFFSET             = 1.2     # perp-offset for extrinsic arrows
LABEL_PADDING           = 0.4     # extra push for extrinsic labels (perp)
PARALLEL_PADDING        = 0.3     # extra push for extrinsic labels (along slanted links)
VERTICAL_LABEL_PAD      = 0.9     # extra push for extrinsic labels on truly vertical links
SCALE_INTRINSIC         = 3       # scaling of intrinsic node layout
INTRINSIC_ARROW_PADDING = 0.8     # how much to trim off each end of intrinsic arrows
CIRCLE_HEAD_RADIUS      = 0.18    # radius of the little circle arrow-heads
CIRCLE_POINTER_SIZE     = CIRCLE_HEAD_RADIUS * 2  # diameter for pointer_size
FIG_WIDTH_PER           = 8       # inches per panel (we’ll have 3 panels per figure)
FIG_HEIGHT              = 8       # inches total height

# **New**: extra trimming for external arrows
EXTRINSIC_TRIM          = 1       # additional data–space units to shorten extrinsic arrows

# You can adjust this dict to offset specific ion channels vertically
ION_CHANNEL_Y_OFFSETS = {
    1: -0.12,  # bring channel 1 down by 0.12
}

# -------------------------------
# Self-loop offsets (u: (x_offset, y_offset))
# -------------------------------
SELF_LOOP_OFFSETS = {
    1: (-0.4, 0.7),  # SS self-loop (default)
    2: (0.7,  0.0),
    3: (0.4,  -0.7),
    4: (0.8, 0.0),   # DP self-loop
}

# -------------------------------
# Self-loop LABEL offsets
# -------------------------------
SELF_LOOP_LABEL_OFFSETS = {
    1: (-0.4, 1.2),  # SS label
    2: (1.0,  0.2),
    3: (0.0,  -1.2),
    4: (1.0,  0.2),
}

# -------------------------------
# Raw mapping (keys will be normalized)
# -------------------------------
RAW_TITLE_MAP = {
    'covariate1': 'baseline',
    'covariate2': 'sustained',
    'covariate3': 'transient'
}

# -------------------------------
# Layout & images
# -------------------------------
region_positions = {
    1: (0, GRID_SPACING),
    2: (GRID_SPACING, GRID_SPACING),
    3: (0, 0),
    4: (GRID_SPACING, 0)
}
intrinsic_base_positions = {
    1: (-1, 0.5),
    2: (0, 2.0),
    3: (1, 1.0),
    4: (0, 0)
}
cell_nodes = {1: "SS", 2: "SP", 3: "II", 4: "DP"}
node_images = {
    1: r"D:\signal_data\fog\pics\Picture6.png",
    2: r"D:\signal_data\fog\pics\Picture5.png",
    3: r"D:\signal_data\fog\pics\Picture7.png",
    4: r"D:\signal_data\fog\pics\Picture8.png"
}
shape_images = {
    1: r"D:\signal_data\fog\pics\Picture10.png",
    2: r"D:\signal_data\fog\pics\Picture20.png",
    3: r"D:\signal_data\fog\pics\Picture30.png"
}

# -------------------------------
# Helper: draw a circular self-loop with a round head
# -------------------------------
def draw_self_loop_with_circle(ax, position, loop_radius=0.4, color=NEG_COLOR,
                               lw=2, offset=(0, 0), arrow_angle_deg=270,
                               pointer_size=CIRCLE_POINTER_SIZE):
    cx, cy = position[0] + offset[0], position[1] + offset[1]
    loop = patches.Circle((cx, cy), loop_radius,
                          edgecolor=color, facecolor='none',
                          lw=lw, zorder=12)
    ax.add_patch(loop)
    angle = np.deg2rad(arrow_angle_deg)
    bx = cx + loop_radius * np.cos(angle)
    by = cy + loop_radius * np.sin(angle)
    head = patches.Circle((bx, by), pointer_size / 2,
                          edgecolor=color, facecolor=color,
                          zorder=13)
    ax.add_patch(head)
    return loop, head

# -------------------------------
# Load .mat & extract Ep/Pnames segments
# -------------------------------
def process_file(mat_file_path):
    import re
    import numpy as np
    import scipy.io as sio
    from scipy.sparse import issparse

    # --- robust regex (spaces allowed; optional 3rd index in T)
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

    # --- load MAT
    mat = sio.loadmat(mat_file_path, squeeze_me=True, struct_as_record=False)
    if 'Results' not in mat:
        raise KeyError(f"'Results' not found in {mat_file_path}")
    res = mat['Results']
    if not hasattr(res, 'PEB_thresholded'):
        raise KeyError(f"'Results.PEB_thresholded' not found in {mat_file_path}")
    peb = res.PEB_thresholded

    # --- flatten Ep / Pnames regardless of save style
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
        Ep_s = Ep[s * Np:(s + 1) * Np] if nseg > 1 else Ep

        cov_ex    = {}  # extrinsic connections per covariate
        cov_int   = {}  # intrinsic (H/… groups) per covariate
        cov_shape = {}  # ion channels (T) per covariate

        for i, pn in enumerate(Pnames):
            val = Ep_s[i]

            # --- ALWAYS capture ion-channel (T) entries, even when val == 0 ---
            m_sh = pat_shape.search(pn)
            if m_sh:
                cov_idx, area, shape_idx = map(int, m_sh.groups())
                key = f"Covariate {cov_idx}"
                cov_shape.setdefault(key, []).append((area, shape_idx, val))
                # Do not 'continue' here; we still might parse extrinsic/intrinsic below if val != 0

            # For extrinsic/intrinsic edges, we keep your original behavior: skip zeros
            if val == 0:
                continue

            # Normalize "Covariate X: ..." split
            if ':' in pn:
                parts = re.split(r'\s*:\s*', pn, maxsplit=1)
                if len(parts) == 2:
                    cov, conn = parts[0].strip(), parts[1].strip()
                else:
                    cov, conn = pn, ''
            else:
                cov, conn = pn, ''

            # 1) EXTRINSIC: A/AN{ins}(tgt,src)
            m_ex = pat_ex.search(conn)
            if m_ex:
                typ, ins, tgt, src = m_ex.groups()
                src_i, tgt_i = int(src), int(tgt)
                label = f"{typ}{ins}"
                cov_ex.setdefault(cov, []).append((src_i, tgt_i, label, val))

            # 2) INTRINSIC: H(...)
            if re.match(r'(?i)^covariate\s+\d+\b', cov) and conn.startswith('H'):
                cov_int.setdefault(cov, []).append({"label": conn, "Ep": val})

        segments.append({
            'cov_extrinsic':    cov_ex,
            'covariate_groups': cov_int,
            'cov_shape':        cov_shape,
            'segment':          s + 1
        })

    return segments


# -------------------------------
# Plot one covariate on one Axes
# -------------------------------
def plot_connectivity_for_covariate(ax, cov, data):
    ax.set_aspect('equal','box')
    ax.axis('off')

    # draw regions
    for pos in region_positions.values():
        ax.add_patch(patches.Circle(pos, REGION_RADIUS,
                                    facecolor='none', edgecolor='black',
                                    lw=0, zorder=1))

    # ion-channel icons & arrows
    cov_shapes = data.get('cov_shape', {}).get(cov, [])
    shapes_by_area = {}
    for area, shape_idx, ep_val in cov_shapes:
        shapes_by_area.setdefault(area, []).append((shape_idx, ep_val))

    base_angles = {1:135, 2:45, 3:225, 4:315}
    spread_deg  = 30
    for area, items in shapes_by_area.items():
        items_sorted = sorted(items, key=lambda x: x[0])
        n = len(items_sorted)
        angles = (np.linspace(base_angles[area]-spread_deg/2,
                              base_angles[area]+spread_deg/2, n)
                  if n>1 else [base_angles[area]])
        if area in (1,3):
            angles = angles[::-1]
        cx, cy = region_positions[area]
        for (shape_idx, ep_val), angle in zip(items_sorted, angles):
            dx = np.cos(np.deg2rad(angle))*(REGION_RADIUS+0.5)
            dy = np.sin(np.deg2rad(angle))*(REGION_RADIUS+0.5)
            x0, y0 = cx+dx, cy+dy
            y0 += ION_CHANNEL_Y_OFFSETS.get(shape_idx, 0)

            try:
                img_arr = plt.imread(shape_images[shape_idx])
                img_box = OffsetImage(img_arr, zoom=SHAPE_ICON_ZOOM)
                ab = AnnotationBbox(img_box, (x0, y0),
                                    frameon=False, zorder=5)
                ax.add_artist(ab)
            except FileNotFoundError:
                ax.scatter(x0, y0, s=100, color='gray', zorder=5)

            if ep_val != 0:
                col = POS_COLOR if ep_val > 0 else NEG_COLOR
                dy_arrow = SHAPE_ARROW_LENGTH if ep_val>0 else -SHAPE_ARROW_LENGTH
                ax.annotate('', xy=(x0, y0+dy_arrow), xytext=(x0, y0),
                            arrowprops=dict(arrowstyle='-|>', color=col, lw=1.5),
                            zorder=6)
                ax.text(x0, y0+dy_arrow, fmt_trunc_str(ep_val, 2),
                        color=col, fontsize=8,
                        ha='center',
                        va='bottom' if ep_val>0 else 'top',
                        bbox=dict(facecolor='white', alpha=0.7),
                        zorder=6)

    # extrinsic edges
    cov_ex = data['cov_extrinsic'].get(cov, [])
    grouped_ex = {}
    for src, tgt, lab, ep in cov_ex:
        grouped_ex.setdefault(tuple(sorted((src, tgt))), []).append((src, tgt, lab, ep))
    for conns in grouped_ex.values():
        for idx, (src, tgt, lab, ep) in enumerate(conns):
            base_idx = idx - (len(conns)-1)/2
            sign     = 1 if src < tgt else -1
            off_fac  = base_idx * sign

            P_src, P_tgt = np.array(region_positions[src]), np.array(region_positions[tgt])
            vec = P_tgt - P_src
            L = np.linalg.norm(vec)
            unit = vec/L if L else np.zeros(2)
            perp = np.array([-unit[1], unit[0]])

            trim = REGION_RADIUS + EXTRINSIC_TRIM
            start = P_src + perp*BASE_OFFSET*off_fac + unit*trim
            end   = P_tgt + perp*BASE_OFFSET*off_fac - unit*trim

            color = POS_COLOR if ep>0 else NEG_COLOR

            ax.annotate('', xy=tuple(end), xytext=tuple(start),
                        arrowprops=dict(arrowstyle='->', lw=3,
                                        color=color, mutation_scale=20),
                        zorder=3)

            mid = (start+end)/2
            off_perp = perp*LABEL_PADDING*off_fac
            pad = VERTICAL_LABEL_PAD if abs(unit[0])<0.1 else PARALLEL_PADDING
            off_para = unit*pad*base_idx if abs(unit[0])<abs(unit[1]) else np.zeros(2)
            lbl_pos = mid + off_perp + off_para
            ax.text(lbl_pos[0], lbl_pos[1],
                    f"{lab} ({fmt_trunc_str(ep, 2)})",
                    color='black', fontsize=10, ha='center', va='center',
                    zorder=4)

    # intrinsic edges + self-loops
    cov_int = data['covariate_groups'].get(cov, [])
    for roi in region_positions:
        G = nx.MultiDiGraph()
        for n in cell_nodes:
            G.add_node(n)
        for ent in cov_int:
            if f",{roi})" in ent['label']:
                d,s,r = map(int, re.findall(r"\((\d+),(\d+),(\d+)\)", ent['label'])[0])
                if r == roi:
                    G.add_edge(s, d, Ep=ent['Ep'])

        cent = np.mean(list(intrinsic_base_positions.values()), axis=0)
        pos = {n: tuple(np.array(region_positions[roi]) +
                        SCALE_INTRINSIC*(np.array(p)-cent))
               for n,p in intrinsic_base_positions.items()}

        nonself, selfe = [], []
        for u,v,d in G.edges(data=True):
            (selfe if u==v else nonself).append((u,v,d))

        grouped = {}
        for u,v,d in nonself:
            grouped.setdefault(tuple(sorted((u,v))), []).append((u,v,d))
        for edges in grouped.values():
            for i,(u,v,d) in enumerate(edges):
                rad = 0.3 + 0.05*(i - (len(edges)-1)/2)
                col = POS_COLOR if d['Ep']>0 else NEG_COLOR
                p0,p1 = np.array(pos[u]), np.array(pos[v])
                vec2, L2 = p1-p0, np.linalg.norm(p1-p0)
                if L2>2*INTRINSIC_ARROW_PADDING:
                    uv2 = vec2/L2
                    p0 += uv2*INTRINSIC_ARROW_PADDING
                    p1 -= uv2*INTRINSIC_ARROW_PADDING
                arrow = FancyArrowPatch(p0, p1,
                                       connectionstyle=f"arc3,rad={rad}",
                                       arrowstyle='-|>' if u!=3 else '-',
                                       mutation_scale=20,
                                       color=col, lw=2, zorder=12)
                ax.add_patch(arrow)
                if u == 3:
                    ax.add_patch(patches.Circle(p1, CIRCLE_HEAD_RADIUS,
                                                facecolor=col, zorder=13))
                path = arrow.get_path().transformed(arrow.get_patch_transform())
                mp = path.interpolated(steps=100).vertices[50]
                ax.text(mp[0], mp[1], fmt_trunc_str(d['Ep'], 2),
                        color=col, fontsize=9, ha='center', va='center',
                        bbox=dict(facecolor='white', alpha=0.7), zorder=13)

        for u,_,d in selfe:
            p = pos[u]
            offset = SELF_LOOP_OFFSETS.get(u, (0,0))
            col = POS_COLOR if d['Ep']>0 else NEG_COLOR
            draw_self_loop_with_circle(ax, p, loop_radius=0.4,
                                       color=col, offset=offset)
            label_offset = SELF_LOOP_LABEL_OFFSETS.get(u, (0,0))
            lbl = (p[0] + label_offset[0], p[1] + label_offset[1])
            ax.text(lbl[0], lbl[1], fmt_trunc_str(d['Ep'], 2),
                    color=col, fontsize=9, ha='center', va='center',
                    bbox=dict(facecolor='white', alpha=0.7), zorder=14)

        for n, coord in pos.items():
            try:
                im = plt.imread(node_images[n])
                icon = OffsetImage(im, zoom=ICON_ZOOM)
                ab   = AnnotationBbox(icon, coord, frameon=False)
                ax.add_artist(ab)
            except FileNotFoundError:
                ax.scatter(*coord, s=200, color='gray', zorder=10)

    xs = [p[0] for p in region_positions.values()]
    ys = [p[1] for p in region_positions.values()]
    pad = REGION_RADIUS + 1
    extra = SHAPE_ARROW_LENGTH + max(abs(v) for v in ION_CHANNEL_Y_OFFSETS.values())
    ax.set_xlim(min(xs) - pad - extra, max(xs) + pad + extra)
    ax.set_ylim(min(ys) - pad - extra, max(ys) + pad + extra)
    ax.margins(0.05)

# -------------------------------
# MAIN entry (variant: 3 figures per segment; each figure = one covariate with iTBS, cTBS, Sham)
# -------------------------------
if __name__ == "__main__":
    # Keep your original paths
    file_paths = [
        r"D:\signal_data\cleaning_git\to_plot\Plot_ctbs_final_BST_n1c.mat",  # cTBS
        r"D:\signal_data\cleaning_git\to_plot\Plot_itbs_final_BST_n1i.mat",  # iTBS
        r"D:\signal_data\cleaning_git\to_plot\Plot_sham_final_BST_n1s.mat",  # Sham
    ]
    all_segments_list = [process_file(fp) for fp in file_paths]

    # Map conditions to segment-lists
    # NOTE: desired order within each figure = iTBS, cTBS, Sham
    condition_segments = {
        'iTBS': all_segments_list[1],
        'cTBS': all_segments_list[0],
        'Sham': all_segments_list[2],
    }
    conditions_order = ['iTBS', 'cTBS', 'Sham']

    # Covariate order & labels
    covariates = [
        ('Covariate 1', 'baseline'),
        ('Covariate 2', 'sustained'),
        ('Covariate 3', 'transient'),
    ]

    outdir = r'D:\signal_data\cleaning_git\to_plot'
    os.makedirs(outdir, exist_ok=True)

    # Assume same number of segments for all; be safe with min length
    num_segments = min(len(v) for v in condition_segments.values())

    for seg_idx in range(num_segments):
        # Build three figures: one per covariate, each with 3 panels (iTBS, cTBS, Sham)
        for cov_key, cov_label in covariates:
            fig, axs = plt.subplots(
                1, len(conditions_order),
                figsize=(FIG_WIDTH_PER * len(conditions_order), FIG_HEIGHT),
                squeeze=False
            )
            axs = axs.flatten()
            fig.patch.set_facecolor('none')
            for ax in axs:
                ax.set_facecolor('none')

            # Plot each condition in requested order
            for j, cond in enumerate(conditions_order):
                seg_data = condition_segments[cond][seg_idx]  # pick the same segment index per condition
                plot_connectivity_for_covariate(axs[j], cov_key, seg_data)
                axs[j].set_title(cond, fontsize=12)

            # Figure title and layout
            fig.suptitle(f"{cov_label.capitalize()} – Segment {seg_idx+1}", fontsize=16)
            fig.subplots_adjust(left=0.05, right=0.98, top=0.9, bottom=0.05, wspace=0.05)

            # Base name for saving
            fname_base = f"{cov_label}_segment{seg_idx+1}"
            svg_path   = os.path.join(outdir, fname_base + '.svg')

            # Save original SVG (transparent like before)
            fig.savefig(svg_path, format='svg', bbox_inches='tight', transparent=True)
            print(f"Saved figure to: {svg_path}")

            # ---- Extra saves: color-mode folders & formats (SVG/TIFF/JPG @ >=300 DPI) ----
            DPI_EXPORT = 600  # ensure >300 DPI

            overtime_root = os.path.join(outdir, 'overtime_small')
            subfolders = ['RBG', 'CMYK', 'Grayscale']  # per your requested names
            for sub in subfolders:
                os.makedirs(os.path.join(overtime_root, sub), exist_ok=True)

            base = fname_base

            # 1) Save SVG copies into each subfolder
            for sub in subfolders:
                svg_copy = os.path.join(overtime_root, sub, base + '.svg')
                fig.savefig(svg_copy, format='svg',
                            bbox_inches='tight', transparent=True, dpi=DPI_EXPORT)

            # 2) High-res raster once in RGB (white background)
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=DPI_EXPORT, bbox_inches='tight', facecolor='white')
            buf.seek(0)
            img_rgb = Image.open(buf).convert('RGB')

            # 3) CMYK and Grayscale versions
            img_cmyk = img_rgb.convert('CMYK')
            img_gray = img_rgb.convert('L')

            # 4) Save TIFF & JPG per color mode

            # RBG
            rbg_dir = os.path.join(overtime_root, 'RBG')
            img_rgb.save(os.path.join(rbg_dir, base + '.tif'), format='TIFF', dpi=(DPI_EXPORT, DPI_EXPORT))
            img_rgb.save(os.path.join(rbg_dir, base + '.jpg'), format='JPEG',
                         quality=95, subsampling=0, optimize=True, dpi=(DPI_EXPORT, DPI_EXPORT))

            # CMYK
            cmyk_dir = os.path.join(overtime_root, 'CMYK')
            img_cmyk.save(os.path.join(cmyk_dir, base + '.tif'), format='TIFF', dpi=(DPI_EXPORT, DPI_EXPORT))
            img_cmyk.save(os.path.join(cmyk_dir, base + '.jpg'), format='JPEG',
                          quality=95, subsampling=0, optimize=True, dpi=(DPI_EXPORT, DPI_EXPORT))

            # Grayscale
            gray_dir = os.path.join(overtime_root, 'Grayscale')
            img_gray.save(os.path.join(gray_dir, base + '.tif'), format='TIFF', dpi=(DPI_EXPORT, DPI_EXPORT))
            img_gray.save(os.path.join(gray_dir, base + '.jpg'), format='JPEG',
                          quality=95, subsampling=0, optimize=True, dpi=(DPI_EXPORT, DPI_EXPORT))
            # -------------------------------------------------------------------------------

            plt.show()
