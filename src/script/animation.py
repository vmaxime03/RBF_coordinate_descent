import math
import os
import re
import glob
import json
import sys
 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import FancyArrowPatch
from matplotlib.widgets import Button
 


# RBF 
def gaussian(r, s):
    return np.exp(-r * r / (2 * s * s))

def gaussian_d(r, s):
    return (-r / (s * s)) * gaussian(r, s)


RBF_MAP = {
    "gaussian": (gaussian, gaussian_d),
}


# SDF
def eval_sdf(X, Y, data):
    rbf_name = data.get("rbf", "gaussian")
    phi, dphi = RBF_MAP[rbf_name]

    result = np.zeros_like(X, dtype=float)
    for p in data["points"]:
        px, py   = p["px"], p["py"]
        alpha    = p["alpha"]
        bx, by   = p["bx"], p["by"]
        s        = p["s"]

        DX = X - px
        DY = Y - py
        n  = np.maximum(np.sqrt(DX**2 + DY**2), 1e-8)

        result += alpha * phi(n, s) + dphi(n, s) * (bx * DX / n + by * DY / n)


    return result


# LOAD FILES
def load_snapshots(directory):
    pattern = os.path.join(directory, "*sdf.json")
    files   = glob.glob(pattern)

    def iteration(path):
        m = re.search(r"(\d+)sdf\.json$", os.path.basename(path))
        return int(m.group(1)) if m else -1

    files = sorted(files, key=iteration)
    snapshots = []
    for f in files:
        with open(f) as fp:
            snapshots.append((iteration(f), json.load(fp)))

    return snapshots



# MAIN
def main():

    snapshots = load_snapshots(sys.argv[1])
    if not snapshots:
        sys.exit(1)

    g   = 2.0;
    res = 300
    x   = np.linspace(-g, g, res)
    y   = np.linspace(-g, g, res)
    X, Y = np.meshgrid(x, y)

    # compute frames
    Zs = [eval_sdf(X, Y, data) for _, data in snapshots]
    abs_max = max(max(abs(Z.min()), abs(Z.max())) for Z in Zs)
 
    # layout
    fig, ax = plt.subplots(figsize=(7, 7))
    fig.patch.set_facecolor("#0f0f1e")
    ax.set_facecolor("#0f0f1e")
    ax.set_aspect("equal")
    ax.set_xlim(-g, g)
    ax.set_ylim(-g, g)
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor("#333")
 
    # initial draw
    def draw_frame(idx):
        ax.cla()
        ax.set_facecolor("#0f0f1e")
        ax.set_aspect("equal")
        ax.set_xlim(-g, g)
        ax.set_ylim(-g, g)
        ax.tick_params(colors="white")
        for spine in ax.spines.values():
            spine.set_edgecolor("#333")
 
        it, data = snapshots[idx]
        Z = Zs[idx]
 
        ax.contourf(X, Y, Z, levels=60,
                    cmap="RdBu", vmin=-abs_max, vmax=abs_max, zorder=1)
        ax.contour(X, Y, Z, levels=20,
                   colors="white", linewidths=0.4, alpha=0.3, zorder=2)
        ax.contour(X, Y, Z, levels=[0],
                   colors=["#ff4444"], linewidths=2.5, zorder=3)


        # polyline TODO read from json
        px, py = [-1, 0, 1], [1, 0, 1]
        ax.plot(px, py, color="lime", linewidth=2, zorder=4)

         
        # control points + beta vectors
        for p in data["points"]:
            ax.scatter(p["px"], p["py"],
                       s=90, c="gold", edgecolors="black", linewidths=0.8, zorder=5)
 
            blen = math.sqrt(p["bx"]**2 + p["by"]**2)
            if blen > 1e-6:
                scale = 0.5 / blen
                arr = FancyArrowPatch(
                    (p["px"], p["py"]),
                    (p["px"] + p["bx"] * scale, p["py"] + p["by"] * scale),
                    arrowstyle="-|>", color="yellow",
                    mutation_scale=14, linewidth=2, zorder=6
                )
                ax.add_patch(arr)
 
        ax.set_title(f"Iteration {it}  [{idx+1}/{len(snapshots)}]",
                     color="white", fontsize=12, pad=8)
 
 
    def update(frame_idx):
        draw_frame(frame_idx)
        return []
 
    ani = animation.FuncAnimation(
        fig, update,
        frames=len(snapshots),
        interval=300,        
        repeat=True,
        blit=False,
    )
 
    if (len(sys.argv) > 2) :
        from matplotlib.animation import PillowWriter

        writer = PillowWriter(fps=3) 
        ani.save(sys.argv[2], writer=writer)
    else : 
        plt.show()


if __name__ == "__main__":
    main()
