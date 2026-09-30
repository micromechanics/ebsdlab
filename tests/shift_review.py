"""Classify pytest-mpl image differences into 'explained by a small shift' and 'real change'.

Usage: python tests/shift_review.py <mpl-results-dir> <output-dir> [maxShift] [tolerance]

See tests/Readme.md, section "Reviewing image differences".

For every test folder containing baseline.png and result.png:
  - try all integer shifts (dy, dx) with |dy|, |dx| <= maxShift
  - for each shift, average the residual over a 7x7 window, so a pixel only counts as
    'explained' if its whole neighbourhood matches under the same shift
  - pixels whose residual under the locally best shift is still > tolerance are 'unexplained'
Writes overlay/shift-map PNGs and index.html. Nothing is approved automatically.
"""
import html
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import uniform_filter  # noqa: E402

resultsDir, outDir = Path(sys.argv[1]), Path(sys.argv[2])
maxShift = int(sys.argv[3]) if len(sys.argv) > 3 else 2
tol = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
outDir.mkdir(parents=True, exist_ok=True)


def analyse(baseline, result):
    """Return per-pixel raw diff, best shift index, residual at best shift and the shift list."""
    shifts = [(dy, dx) for dy in range(-maxShift, maxShift + 1) for dx in range(-maxShift, maxShift + 1)]
    resid = np.empty((len(shifts),) + baseline.shape[:2], dtype=np.float32)
    smooth = np.empty_like(resid)
    for i, (dy, dx) in enumerate(shifts):
        resid[i] = np.abs(baseline - np.roll(result, (dy, dx), axis=(0, 1))).max(-1)
        smooth[i] = uniform_filter(resid[i], 7)
    best = smooth.argmin(0)
    bestResid = np.take_along_axis(resid, best[None], 0)[0]
    raw = resid[shifts.index((0, 0))]
    return raw, best, bestResid, shifts


rows = []
for folder in sorted(p for p in resultsDir.iterdir() if p.is_dir()):
    basePath, resPath = folder / 'baseline.png', folder / 'result.png'
    name = folder.name
    if not resPath.exists():
        continue
    target = outDir / name
    target.mkdir(exist_ok=True)
    result = mpimg.imread(resPath)[..., :3]
    plt.imsave(target / 'result.png', result)
    if not basePath.exists():
        rows.append({'name': name, 'status': 'no baseline'})
        continue
    baseline = mpimg.imread(basePath)[..., :3]
    plt.imsave(target / 'baseline.png', baseline)
    if baseline.shape != result.shape:
        rows.append({'name': name, 'status': f'size differs {baseline.shape} vs {result.shape}'})
        continue
    raw, best, bestResid, shifts = analyse(baseline, result)
    changed = raw > tol
    unexplained = bestResid > tol
    shifted = changed & ~unexplained

    # overlay: dimmed baseline, green = explained by shift, red = unexplained
    gray = baseline.mean(-1, keepdims=True) * 0.35 + 0.65
    overlay = np.repeat(gray, 3, -1)
    overlay[shifted] = [0.1, 0.75, 0.2]
    overlay[unexplained] = [0.9, 0.1, 0.1]
    # zoom window: the zoomW x zoomH region with the most unexplained (red) pixels;
    # fall back to the most changed pixels if nothing is unexplained
    zoomW, zoomH, zoom = 160, 120, 5
    focus = unexplained if unexplained.any() else changed
    density = uniform_filter(focus.astype(np.float32), (zoomH, zoomW))
    cy, cx = np.unravel_index(density.argmax(), density.shape)
    y0 = int(np.clip(cy - zoomH // 2, 0, raw.shape[0] - zoomH))
    x0 = int(np.clip(cx - zoomW // 2, 0, raw.shape[1] - zoomW))
    for label, img in (('baseline', baseline), ('result', result), ('overlay', overlay)):
        crop = np.clip(img[y0:y0 + zoomH, x0:x0 + zoomW], 0, 1)
        plt.imsave(target / f'zoom_{label}.png', np.repeat(np.repeat(crop, zoom, 0), zoom, 1))
    zoomBox = {'x': x0, 'y': y0, 'w': zoomW, 'h': zoomH, 'zoom': zoom,
                          'what': 'unexplained (red) pixels' if unexplained.any() else 'changed pixels'}
    boxed = overlay.copy()
    boxed[[y0, y0 + zoomH - 1], x0:x0 + zoomW] = [0.1, 0.3, 0.9]
    boxed[y0:y0 + zoomH, [x0, x0 + zoomW - 1]] = [0.1, 0.3, 0.9]
    plt.imsave(target / 'overlay.png', np.clip(boxed, 0, 1))
    plt.imsave(target / 'diff.png', np.clip(raw * 5, 0, 1), cmap='gray')

    # shift map: dominant shift per 25px tile, only where something changed
    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.imshow(gray[..., 0], cmap='gray', vmin=0, vmax=1)
    tile = 25
    for y0 in range(0, raw.shape[0], tile):
        for x0 in range(0, raw.shape[1], tile):
            m = shifted[y0:y0 + tile, x0:x0 + tile]
            if m.sum() < 5:
                continue
            idx = np.bincount(best[y0:y0 + tile, x0:x0 + tile][m]).argmax()
            dy, dx = shifts[idx]
            # arrow points from baseline position to where the result content sits
            ax.arrow(x0 + tile / 2, y0 + tile / 2, -dx * 5, -dy * 5, color='tab:blue', width=1.2,
                              head_width=5, length_includes_head=True)
    ax.set_axis_off()
    ax.set_title('dominant shift per tile (arrow length x5)')
    fig.savefig(target / 'shiftmap.png', bbox_inches='tight')
    plt.close(fig)

    hist = {}
    for i in np.unique(best[shifted]):
        hist[str(shifts[i])] = int((best[shifted] == i).sum())
    nPix = raw.size
    rows.append({'name': name, 'status': 'ok', 'changed': int(changed.sum()),
                              'shifted': int(shifted.sum()), 'unexplained': int(unexplained.sum()),
                              'unexplainedPct': 100 * unexplained.sum() / nPix, 'zoom': zoomBox,
                              'hist': dict(sorted(hist.items(), key=lambda kv: -kv[1])[:5])})

(outDir / 'summary.json').write_text(json.dumps(rows, indent=1))

cards = []
for r in rows:
    n = html.escape(r['name'])
    if r['status'] != 'ok':
        cards.append(f'<section><h2>{n}</h2><p class="warn">{html.escape(r["status"])} &mdash; review '
                                  f'result by eye.</p><img src="{n}/result.png"></section>')
        continue
    verdict = ('likely shift-only' if r['unexplainedPct'] < 0.05 else
                          'contains real changes' if r['unexplainedPct'] > 0.5 else 'minor unexplained changes')
    cls = {'likely shift-only': 'good', 'contains real changes': 'bad'}.get(verdict, 'mid')
    hist = ', '.join(f'{k}: {v}' for k, v in r['hist'].items()) or '&ndash;'
    z = r['zoom']
    cards.append(f'''<section><h2>{n} <span class="badge {cls}">{verdict}</span></h2>
<p>changed px: {r["changed"]} &middot; <b class="g">explained by &le;{maxShift}px shift: {r["shifted"]}</b>
  &middot; <b class="r">unexplained: {r["unexplained"]} ({r["unexplainedPct"]:.3f}%)</b><br>
most common shifts (dy, dx) px: {hist}</p>
<div class="grid">
  <figure><div class="blink" data-a="{n}/baseline.png" data-b="{n}/result.png"><img src="{n}/baseline.png"></div>
    <figcaption>click to toggle baseline/result &middot; <label><input type="checkbox" class="auto"> auto-blink</label>
    &middot; showing: <span class="which">baseline</span></figcaption></figure>
  <figure><img src="{n}/overlay.png"><figcaption>green = shift, red = unexplained, blue box = zoom area</figcaption></figure>
  <figure><img src="{n}/shiftmap.png"><figcaption>shift direction per tile</figcaption></figure>
  <figure><img src="{n}/diff.png"><figcaption>raw |baseline&minus;result| &times;5</figcaption></figure>
</div>
<h3>Zoom &times;{z["zoom"]} on area with most {z["what"]} (x {z["x"]}&ndash;{z["x"] + z["w"]}, y {z["y"]}&ndash;{z["y"] + z["h"]} px)</h3>
<div class="grid3">
  <figure><img src="{n}/zoom_baseline.png"><figcaption>before (baseline)</figcaption></figure>
  <figure><img src="{n}/zoom_overlay.png"><figcaption>green = shift, red = unexplained</figcaption></figure>
  <figure><img src="{n}/zoom_result.png"><figcaption>after (result)</figcaption></figure>
</div></section>''')

(outDir / 'index.html').write_text(f'''<!doctype html><html><head><meta charset="utf-8">
<title>Shift review</title><style>
body{{font-family:sans-serif;margin:16px;background:#fafafa;color:#222}}
section{{background:#fff;border:1px solid #ddd;border-radius:6px;padding:12px;margin-bottom:20px}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(380px,1fr));gap:10px}}
.grid3{{display:grid;grid-template-columns:repeat(3,1fr);gap:10px}}
@media (max-width:900px){{.grid3{{grid-template-columns:1fr}}}}
figure{{margin:0}} img{{width:100%;border:1px solid #ccc;image-rendering:pixelated}}
.blink{{cursor:pointer}} .badge{{font-size:.7em;padding:2px 8px;border-radius:10px;color:#fff}}
.good{{background:#2a8a3a}} .bad{{background:#c0392b}} .mid{{background:#d68910}}
.g{{color:#2a8a3a}} .r{{color:#c0392b}} .warn{{color:#c0392b}}
</style></head><body>
<h1>Shift review (max shift {maxShift}px, tolerance {tol})</h1>
<p>Green pixels differ from the baseline but match it again after a consistent local shift of at
most {maxShift}px. Red pixels match under no such shift: look there for real changes.
Zoom in the browser (Ctrl +) for pixel detail. This page approves nothing.</p>
{"".join(cards)}
<script>
document.querySelectorAll('.blink').forEach(d=>{{
    const img=d.querySelector('img'), fig=d.parentElement, which=fig.querySelector('.which');
    let b=false, timer=null;
    const flip=()=>{{b=!b; img.src=b?d.dataset.b:d.dataset.a; which.textContent=b?'result':'baseline';}};
    d.onclick=flip;
    fig.querySelector('.auto').onchange=e=>{{clearInterval(timer); if(e.target.checked) timer=setInterval(flip,600);}};
}});
</script></body></html>''')
for r in rows:
    print(r['name'], r['status'], f"unexplained={r.get('unexplainedPct', float('nan')):.3f}%",
                f"shifted={r.get('shifted', '-')}", r.get('hist', ''))
