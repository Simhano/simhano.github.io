/* main.js — page interactions: Figure 1 (live beam), Figure 2 (configuration
   diagram), Figure 3 (inverse-unloading demo), nav.
   Depends on fem.js. */
(function () {
  'use strict';
  const { meshDomain, sampleCircle, sampleSegment, Body, banded } = window.FEM;
  const BANDS = 12;
  const INK = '#1C1B19', GRAPHITE = '#5F5C55', VERMILION = '#B83A22', BLUEPRINT = '#2E4A7A', RULE = '#E2DED3';
  const MONO = '"IBM Plex Mono", ui-monospace, monospace';
  const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;

  // ---------------------------------------------------------------------------
  // Canvas helpers
  // ---------------------------------------------------------------------------

  /** Keeps a canvas sized to its wrapper (DPR aware) and maps model -> pixels. */
  function makeView(canvas, box) {
    const ctx = canvas.getContext('2d');
    const view = { ctx, canvas, box, s: 1, ox: 0, oy: 0, w: 0, h: 0 };
    view.resize = () => {
      const r = canvas.parentElement.getBoundingClientRect();
      const dpr = Math.min(window.devicePixelRatio || 1, 2);
      canvas.width = Math.round(r.width * dpr);
      canvas.height = Math.round(r.height * dpr);
      view.w = r.width; view.h = r.height;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const [x0, y0, x1, y1] = view.box;
      view.s = Math.min(r.width / (x1 - x0), r.height / (y1 - y0));
      view.ox = (r.width - view.s * (x0 + x1)) / 2;
      view.oy = (r.height + view.s * (y0 + y1)) / 2;
    };
    view.px = (x) => view.ox + view.s * x;
    view.py = (y) => view.oy - view.s * y;
    view.toModel = (cx, cy) => [(cx - view.ox) / view.s, (view.oy - cy) / view.s];
    new ResizeObserver(() => { view.resize(); view.onResize && view.onResize(); }).observe(canvas.parentElement);
    view.resize();
    return view;
  }

  /** Clip a polygon [[x, y, v], ...] to the half-space sign*(v - level) >= 0. */
  function clip(poly, level, sign) {
    const out = [];
    for (let q = 0; q < poly.length; q++) {
      const A = poly[q], B = poly[(q + 1) % poly.length];
      const da = sign * (A[2] - level), db = sign * (B[2] - level);
      if (da >= 0) out.push(A);
      if ((da >= 0) !== (db >= 0)) {
        const t = da / (da - db);
        out.push([A[0] + t * (B[0] - A[0]), A[1] + t * (B[1] - A[1]), level]);
      }
    }
    return out;
  }

  /**
   * Banded contour plot of nodal-averaged von Mises stress: each triangle is
   * split into iso-bands of the linearly interpolated field, one path per band.
   */
  function drawContours(view, body, x, vmax, { wire = true, alpha = 1 } = {}) {
    const { ctx } = view;
    const { T, m } = body.mesh;
    const vn = body._vn = body.nodal(body.mises, body._vn);
    const paths = Array.from({ length: BANDS }, () => new Path2D());
    const band = (v) => Math.max(0, Math.min(BANDS - 1, Math.floor((v / vmax) * BANDS)));
    for (let e = 0; e < m; e++) {
      const tri = [T[3 * e], T[3 * e + 1], T[3 * e + 2]].map((i) => [view.px(x[2 * i]), view.py(x[2 * i + 1]), vn[i]]);
      const k0 = band(Math.min(tri[0][2], tri[1][2], tri[2][2]));
      const k1 = band(Math.max(tri[0][2], tri[1][2], tri[2][2]));
      for (let k = k0; k <= k1; k++) {
        let poly = tri;
        if (k > 0) poly = clip(poly, (k / BANDS) * vmax, 1);
        if (k < BANDS - 1) poly = clip(poly, ((k + 1) / BANDS) * vmax, -1);
        if (poly.length < 3) continue;
        const p = paths[k];
        p.moveTo(poly[0][0], poly[0][1]);
        for (let q = 1; q < poly.length; q++) p.lineTo(poly[q][0], poly[q][1]);
        p.closePath();
      }
    }
    ctx.globalAlpha = alpha;
    paths.forEach((p, k) => {
      ctx.fillStyle = banded((k + 0.5) / BANDS, BANDS);
      ctx.fill(p);
      ctx.strokeStyle = ctx.fillStyle; // hide hairline seams between triangles
      ctx.lineWidth = 0.6;
      ctx.stroke(p);
    });
    ctx.globalAlpha = 1;
    if (wire) drawEdges(view, body.mesh.edges, x, 'rgba(28,27,25,0.32)', 0.55);
    drawEdges(view, body.mesh.boundaryEdges, x, INK, 1.3);
  }

  function drawEdges(view, E, x, color, width, dash) {
    const { ctx } = view;
    ctx.beginPath();
    for (let q = 0; q < E.length; q += 2) {
      const a = E[q], b = E[q + 1];
      ctx.moveTo(view.px(x[2 * a]), view.py(x[2 * a + 1]));
      ctx.lineTo(view.px(x[2 * b]), view.py(x[2 * b + 1]));
    }
    ctx.strokeStyle = color;
    ctx.lineWidth = width;
    ctx.setLineDash(dash || []);
    ctx.stroke();
    ctx.setLineDash([]);
  }

  /** Classic support symbol: a hatched wall along x = const. */
  function drawWall(view, x, y0, y1) {
    const { ctx } = view;
    const X = view.px(x), Y0 = view.py(y1), Y1 = view.py(y0);
    ctx.strokeStyle = INK;
    ctx.lineWidth = 1.6;
    ctx.beginPath(); ctx.moveTo(X, Y0); ctx.lineTo(X, Y1); ctx.stroke();
    ctx.lineWidth = 0.9;
    ctx.beginPath();
    for (let y = Y0; y <= Y1 - 8; y += 8) { ctx.moveTo(X, y); ctx.lineTo(X - 8, y + 8); }
    ctx.stroke();
  }

  function drawArrow(ctx, x0, y0, x1, y1, color, head = 8) {
    const a = Math.atan2(y1 - y0, x1 - x0);
    ctx.strokeStyle = ctx.fillStyle = color;
    ctx.lineWidth = 1.6;
    ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x1 - head * Math.cos(a - 0.38), y1 - head * Math.sin(a - 0.38));
    ctx.lineTo(x1 - head * Math.cos(a + 0.38), y1 - head * Math.sin(a + 0.38));
    ctx.closePath(); ctx.fill();
  }

  function drawGravity(view, on, x, y) {
    const { ctx } = view;
    ctx.globalAlpha = on ? 1 : 0.3;
    drawArrow(ctx, x, y, x, y + 30, INK, 7);
    ctx.font = `italic 15px "Source Serif 4", Georgia, serif`;
    ctx.fillText('g', x + 7, y + 18);
    ctx.globalAlpha = 1;
  }

  function buildLegend(el, vmax, title) {
    let cells = '';
    for (let k = 0; k < BANDS; k++) cells += `<i style="background:${banded((k + 0.5) / BANDS)}"></i>`;
    el.innerHTML = `<span class="cb-title">${title}</span><span class="cb-bar">${cells}</span>` +
      `<span class="cb-ticks"><span>0</span><span>${(vmax / 2).toExponential(1)}</span><span>${vmax.toExponential(1)}</span></span>`;
  }

  /** Runs `tick` with requestAnimationFrame only while `el` is on screen. */
  function whileVisible(el, tick) {
    let raf = 0, visible = false;
    const loop = () => { tick(); raf = visible ? requestAnimationFrame(loop) : 0; };
    new IntersectionObserver(([entry]) => {
      visible = entry.isIntersecting;
      if (visible && !raf) raf = requestAnimationFrame(loop);
    }).observe(el);
  }

  // ---------------------------------------------------------------------------
  // Figure 1: cantilever with transverse holes under gravity, draggable
  // ---------------------------------------------------------------------------

  function initHero() {
    const canvas = document.getElementById('hero-canvas');
    if (!canvas) return;
    const L = 3, H = 1, h = 0.085;
    const holes = [[0.85, 0.5, 0.2], [1.6, 0.5, 0.2], [2.35, 0.5, 0.2]];
    const boundary = [
      ...sampleSegment(0, 0, L, 0, h), ...sampleSegment(L, 0, L, H, h),
      ...sampleSegment(L, H, 0, H, h), ...sampleSegment(0, H, 0, 0, h),
      ...holes.flatMap(([cx, cy, r]) => sampleCircle(cx, cy, r, h * 0.9)),
    ];
    const inside = (x, y) => x >= -1e-9 && x <= L + 1e-9 && y >= -1e-9 && y <= H + 1e-9 &&
      holes.every(([cx, cy, r]) => (x - cx) ** 2 + (y - cy) ** 2 > r * r);
    const mesh = meshDomain({ boundary, inside, bbox: [0, 0, L, H], h });
    const body = new Body(mesh, { mu: 1, lambda: 4, rho: 1, eta: 0.03 });
    body.fix((x) => x < 1e-6);

    const view = makeView(canvas, [-0.35, -0.85, L + 0.4, H + 0.4]);
    const legend = document.getElementById('hero-legend');
    const status = document.getElementById('hero-status');
    const hint = document.getElementById('hero-hint');
    const ui = {
      gravity: document.getElementById('hero-gravity'),
      stiff: document.getElementById('hero-stiff'),
      stiffVal: document.getElementById('hero-stiff-val'),
      wire: document.getElementById('hero-wire'),
      ghost: document.getElementById('hero-ghost'),
      reset: document.getElementById('hero-reset'),
    };
    const G0 = 0.0065;
    let gravityOn = true, stiffness = 1, t = 0, vmax = 0.05, legendMax = -1;
    let grab = null; // { nodes: [[i, w], ...], target: [x, y], anchor: i }

    if (reduceMotion) body.relax(0, -G0, { maxSteps: 30000, tol: 1e-7 });

    ui.gravity.addEventListener('click', () => {
      gravityOn = !gravityOn;
      ui.gravity.classList.toggle('is-on', gravityOn);
      ui.gravity.setAttribute('aria-pressed', String(gravityOn));
      ui.gravity.textContent = gravityOn ? 'Gravity on' : 'Gravity off';
    });
    ui.stiff.addEventListener('input', () => {
      stiffness = Math.pow(10, +ui.stiff.value);
      body.mu = stiffness; body.lambda = 4 * stiffness;
      ui.stiffVal.textContent = `E × ${stiffness.toFixed(2)}`;
    });
    ui.reset.addEventListener('click', () => { body.reset(); t = 0; });

    const pointer = (ev) => {
      const r = canvas.getBoundingClientRect();
      return view.toModel(ev.clientX - r.left, ev.clientY - r.top);
    };
    canvas.addEventListener('pointerdown', (ev) => {
      ev.preventDefault(); // keep a drag from selecting caption text
      const [mx, my] = pointer(ev);
      const i = body.nearestNode(mx, my, 0.3);
      if (i < 0 || body.fixed[i]) return;
      const R = 0.16, nodes = [];
      for (let k = 0; k < mesh.n; k++) {
        if (body.fixed[k]) continue;
        const d = Math.hypot(body.x[2 * k] - body.x[2 * i], body.x[2 * k + 1] - body.x[2 * i + 1]);
        if (d < R) nodes.push([k, 1 - d / R]);
      }
      grab = { nodes, target: [mx, my], anchor: i };
      canvas.setPointerCapture(ev.pointerId);
      canvas.classList.add('dragging');
      hint.style.opacity = 0;
    });
    canvas.addEventListener('pointermove', (ev) => { if (grab) grab.target = pointer(ev); });
    const release = () => { grab = null; canvas.classList.remove('dragging'); };
    canvas.addEventListener('pointerup', release);
    canvas.addEventListener('pointercancel', release);

    const dragForce = (f) => {
      if (!grab) return;
      const a = grab.anchor, K = 0.12 * stiffness;
      const dx = grab.target[0] - body.x[2 * a], dy = grab.target[1] - body.x[2 * a + 1];
      let wsum = 0;
      for (const [, w] of grab.nodes) wsum += w;
      for (const [k, w] of grab.nodes) {
        f[2 * k] += K * dx * w / wsum;
        f[2 * k + 1] += K * dy * w / wsum;
      }
    };

    function frame() {
      const dt = body.stableDt(0.45);
      const sub = Math.min(160, Math.ceil(0.5 / dt));
      const gy = gravityOn ? -G0 : 0;
      for (let s = 0; s < sub; s++) {
        body.step(dt, { gy, damping: 0.06, extra: dragForce, withStress: s === sub - 1 });
      }
      t += sub * dt;

      let peak = 0;
      for (let e = 0; e < mesh.m; e++) if (body.mises[e] > peak) peak = body.mises[e];
      const targetMax = Math.max(0.02 * stiffness, peak);
      vmax += (targetMax - vmax) * 0.08;
      if (Math.abs(vmax - legendMax) / vmax > 0.04) { legendMax = vmax; buildLegend(legend, vmax, 'Mises stress'); }

      const { ctx } = view;
      ctx.clearRect(0, 0, view.w, view.h);
      drawWall(view, 0, -0.08, H + 0.08);
      drawContours(view, body, body.x, vmax, { wire: ui.wire.checked });
      if (ui.ghost.checked) {
        drawEdges(view, mesh.boundaryEdges, mesh.X, INK, 1.2, [5, 4]);
        ctx.font = `italic 14px "Source Serif 4", Georgia, serif`;
        ctx.fillStyle = INK;
        ctx.textAlign = 'right';
        ctx.fillText('unloaded shape', view.px(L), view.py(H) - 8);
        ctx.textAlign = 'left';
      }
      drawGravity(view, gravityOn, 26, 16);
      ctx.font = `14px "Source Serif 4", Georgia, serif`;
      ctx.fillStyle = GRAPHITE;
      ctx.fillText(gravityOn ? 'on Earth' : 'gravity off', 50, 30);
      if (grab) {
        const a = grab.anchor;
        const x0 = view.px(body.x[2 * a]), y0 = view.py(body.x[2 * a + 1]);
        const x1 = view.px(grab.target[0]), y1 = view.py(grab.target[1]);
        if (Math.hypot(x1 - x0, y1 - y0) > 12) drawArrow(ctx, x0, y0, x1, y1, VERMILION);
        ctx.fillStyle = VERMILION;
        ctx.beginPath(); ctx.arc(x0, y0, 3.5, 0, 7); ctx.fill();
      }
      let umax = 0;
      for (let i = 0; i < mesh.n; i++) {
        umax = Math.max(umax, Math.hypot(body.x[2 * i] - mesh.X[2 * i], body.x[2 * i + 1] - mesh.X[2 * i + 1]));
      }
      status.textContent = `${mesh.n} nodes · ${mesh.m} linear triangles · t = ${t.toFixed(0)} · max |u| = ${umax.toFixed(3)}`;
    }
    whileVisible(canvas, frame);
  }

  // ---------------------------------------------------------------------------
  // Figure 3: recover the stress-free shape of a gravity-loaded body
  // ---------------------------------------------------------------------------

  /** Fill every triangle of a mesh, placed at positions x, with one color. */
  function fillMesh(view, mesh, x, color) {
    const { T, m } = mesh, p = new Path2D();
    for (let e = 0; e < m; e++) {
      const i = T[3 * e], j = T[3 * e + 1], k = T[3 * e + 2];
      p.moveTo(view.px(x[2 * i]), view.py(x[2 * i + 1]));
      p.lineTo(view.px(x[2 * j]), view.py(x[2 * j + 1]));
      p.lineTo(view.px(x[2 * k]), view.py(x[2 * k + 1]));
      p.closePath();
    }
    view.ctx.fillStyle = color;
    view.ctx.fill(p);
  }

  function initDemo() {
    const canvasA = document.getElementById('demo-canvas-a');
    const canvasB = document.getElementById('demo-canvas-b');
    if (!canvasA || !canvasB) return;
    const a = 1.0, b = 1.1, h = 0.09, G = 0.8, TOL = 1e-5, MAX_IT = 80, DONE = 2e-4;
    const bnd = [];
    const K = Math.round(Math.PI * Math.sqrt((a * a + b * b) / 2) / h);
    for (let i = 0; i <= K; i++) {
      const t = -Math.PI / 2 + Math.PI * i / K;
      bnd.push([a * Math.cos(t), b * Math.sin(t)]);
    }
    bnd.push(...sampleSegment(0, b, 0, -b, h).slice(1));
    const mesh = meshDomain({
      boundary: bnd, bbox: [0, -b, a, b], h,
      inside: (x, y) => x >= -1e-9 && (x / a) ** 2 + (y / b) ** 2 <= 1 + 1e-9,
    });
    const body = new Body(mesh, { mu: 1, lambda: 10, rho: 1, eta: 0 });
    body.fix((x) => x < 1e-6);

    const box = [-0.3, -1.75, 1.3, 1.25];
    const viewA = makeView(canvasA, box), viewB = makeView(canvasB, box);
    // Every other free boundary node carries a mismatch / update arrow.
    const arrowNodes = [...new Set(mesh.boundaryEdges)].filter((i) => !body.fixed[i]).filter((_, q) => q % 2 === 0);
    const el = {
      run: document.getElementById('demo-run'),
      step: document.getElementById('demo-step'),
      reset: document.getElementById('demo-reset'),
      truth: document.getElementById('demo-truth'),
      alpha: document.getElementById('demo-alpha'),
      alphaVal: document.getElementById('demo-alpha-val'),
      text: document.getElementById('demo-text'),
      status: document.getElementById('demo-status'),
      chart: document.getElementById('conv-chart'),
      legend: document.getElementById('demo-legend'),
    };
    const height = 2 * b;
    const mis = new Float64Array(2 * mesh.n); // x_obs - x_pred from the last solve
    let xObs = null, Xu = null, phase = 'observe', hasMis = false, auto = false, budget = 0;
    let history = []; // [{ res, geo }] per solve; index = number of updates applied
    let vmax = 2.3, alpha = +el.alpha.value; // vmax is reset from the scan's peak stress

    const say = (html) => { el.text.innerHTML = html; };
    const setButtons = () => {
      const busy = phase === 'observe' || phase === 'solving';
      el.run.disabled = (busy && !auto) || phase === 'done' || phase === 'observe';
      el.step.disabled = busy || phase === 'done';
      el.run.textContent = auto && phase !== 'done' ? 'Pause' : 'Solve';
    };

    function start() {
      body.setReference(mesh.X);
      body.reset();
      xObs = null; Xu = null; history = []; hasMis = false; auto = false;
      mis.fill(0);
      phase = 'observe';
      say('Generating a synthetic scan: placing the true stress-free shape under gravity…');
      setButtons();
      drawChart();
    }

    function beginIteration() {
      // Move the guess by a fraction of the last mismatch: X <- X + alpha (x_obs - x_pred)
      if (hasMis) {
        for (let i = 0; i < mesh.n; i++) {
          if (body.fixed[i]) continue;
          Xu[2 * i] += alpha * mis[2 * i];
          Xu[2 * i + 1] += alpha * mis[2 * i + 1];
        }
      }
      // Re-solve with the new reference; warm start from the last solution.
      const prev = body.Xref;
      body.setReference(Xu);
      for (let i = 0; i < body.x.length; i++) body.x[i] += Xu[i] - prev[i];
      body.v.fill(0);
      phase = 'solving';
      setButtons();
    }

    function finishIteration() {
      let res = 0, geo = 0;
      for (let i = 0; i < mesh.n; i++) {
        if (body.fixed[i]) { mis[2 * i] = mis[2 * i + 1] = 0; continue; }
        mis[2 * i] = xObs[2 * i] - body.x[2 * i];
        mis[2 * i + 1] = xObs[2 * i + 1] - body.x[2 * i + 1];
        res = Math.max(res, Math.hypot(mis[2 * i], mis[2 * i + 1]));
        geo = Math.max(geo, Math.hypot(Xu[2 * i] - mesh.X[2 * i], Xu[2 * i + 1] - mesh.X[2 * i + 1]));
      }
      hasMis = true;
      history.push({ res: res / height, geo: geo / height });
      const k = history.length - 1;
      const pct = (100 * res / height).toFixed(res / height < 0.01 ? 3 : 1);
      if (k === 0) {
        say(`<b>The naive guess.</b> Treat the scan as if it were stress-free and apply gravity again. ` +
            `The body sags a second time; the mismatch in (b) is ${pct}% of its height. ` +
            `The blue arrows in (a) show the correction: α = ${alpha.toFixed(2)} of that mismatch. ` +
            `Press <b>Solve</b> or <b>Step</b>.`);
        phase = 'ready';
      } else if (res / height < DONE || k >= MAX_IT) {
        say(`<b>Converged after ${k} updates.</b> Under gravity, the guess in (a) now sags into the scanned shape in (b) ` +
            `to within ${pct}% of its height. Show the true stress-free shape to compare, or reset and try a larger α.`);
        phase = 'done'; auto = false;
      } else {
        say(`<b>Update ${k}.</b> The guess in (a) moved by α × mismatch, and gravity was applied again. ` +
            `Remaining mismatch in (b): ${pct}% of the height.`);
        phase = 'ready';
      }
      setButtons();
      drawChart();
      if (auto) budget = 6; // short pause so each update is visible
    }

    el.run.addEventListener('click', () => {
      if (auto) { auto = false; setButtons(); return; }
      auto = true;
      if (phase === 'ready') beginIteration();
      setButtons();
    });
    el.step.addEventListener('click', () => { auto = false; if (phase === 'ready') beginIteration(); });
    el.reset.addEventListener('click', start);
    el.truth.addEventListener('change', draw);
    el.alpha.addEventListener('input', () => {
      alpha = +el.alpha.value;
      el.alphaVal.textContent = alpha.toFixed(2);
      draw();
    });

    function drawChart() {
      const c = el.chart, ctx = c.getContext('2d');
      const r = c.getBoundingClientRect(), dpr = Math.min(window.devicePixelRatio || 1, 2);
      c.width = r.width * dpr; c.height = r.height * dpr;
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      const W = r.width, Hh = r.height, padL = 40, padB = 18, padT = 6;
      ctx.clearRect(0, 0, W, Hh);
      const lo = -5, hi = 0; // log10 range
      const xMax = Math.max(10, Math.ceil((history.length - 1) / 10) * 10);
      const X = (k) => padL + (W - padL - 8) * k / xMax;
      const Y = (v) => padT + (Hh - padT - padB) * (hi - Math.max(lo, Math.log10(v))) / (hi - lo);
      ctx.font = `11px ${MONO}`;
      ctx.fillStyle = GRAPHITE;
      ctx.strokeStyle = RULE;
      ctx.lineWidth = 1;
      for (let p = lo; p <= hi; p++) {
        const y = Y(10 ** p);
        ctx.beginPath(); ctx.moveTo(padL, y); ctx.lineTo(W - 8, y); ctx.stroke();
        ctx.fillText(p === 0 ? '1' : `1e${p}`, 0, y + 4);
      }
      const tick = xMax <= 20 ? 5 : 10;
      for (let k = 0; k <= xMax; k += tick) ctx.fillText(String(k), X(k) - 3, Hh - 3);
      const line = (key, color, dash) => {
        if (!history.length) return;
        ctx.strokeStyle = ctx.fillStyle = color; ctx.lineWidth = 1.6; ctx.setLineDash(dash);
        ctx.beginPath();
        history.forEach((hh, k) => (k ? ctx.lineTo : ctx.moveTo).call(ctx, X(k), Y(hh[key])));
        ctx.stroke(); ctx.setLineDash([]);
        history.forEach((hh, k) => { ctx.beginPath(); ctx.arc(X(k), Y(hh[key]), 2.2, 0, 7); ctx.fill(); });
      };
      line('geo', BLUEPRINT, [2, 3]);
      line('res', VERMILION, [5, 3]);
    }

    /** Arrows from points p (+ scale * d) at the arrow nodes, skipping tiny ones. */
    function drawVectors(view, p, d, scale, color) {
      for (const i of arrowNodes) {
        const x0 = view.px(p[2 * i]), y0 = view.py(p[2 * i + 1]);
        const x1 = view.px(p[2 * i] + scale * d[2 * i]), y1 = view.py(p[2 * i + 1] + scale * d[2 * i + 1]);
        if (Math.hypot(x1 - x0, y1 - y0) > 4) drawArrow(view.ctx, x0, y0, x1, y1, color, 5);
      }
    }

    function draw() {
      const showArrows = hasMis && phase !== 'solving' && phase !== 'observe';

      // (a) the stress-free guess, which is what we optimize
      let ctx = viewA.ctx;
      ctx.clearRect(0, 0, viewA.w, viewA.h);
      drawWall(viewA, 0, -b - 0.1, b + 0.1);
      if (Xu) {
        fillMesh(viewA, mesh, Xu, 'rgba(46,74,122,0.08)');
        drawEdges(viewA, mesh.edges, Xu, 'rgba(46,74,122,0.25)', 0.6);
        drawEdges(viewA, mesh.boundaryEdges, Xu, BLUEPRINT, 1.8);
        if (showArrows && phase !== 'done') drawVectors(viewA, Xu, mis, alpha, BLUEPRINT);
      }
      if (el.truth.checked) drawEdges(viewA, mesh.boundaryEdges, mesh.X, INK, 1.6, [2, 4]);
      ctx.font = `italic 14px "Source Serif 4", Georgia, serif`;
      ctx.fillStyle = GRAPHITE;
      ctx.textAlign = 'right';
      ctx.fillText('no gravity', viewA.w - 12, 24);
      ctx.textAlign = 'left';

      // (b) the guess under gravity, compared with the scan
      ctx = viewB.ctx;
      ctx.clearRect(0, 0, viewB.w, viewB.h);
      drawWall(viewB, 0, -b - 0.1, b + 0.1);
      drawContours(viewB, body, body.x, vmax, { wire: true, alpha: phase === 'observe' ? 0.5 : 0.9 });
      if (xObs) drawEdges(viewB, mesh.boundaryEdges, xObs, VERMILION, 2, [7, 5]);
      if (showArrows) drawVectors(viewB, body.x, mis, 1, VERMILION);
      drawGravity(viewB, true, viewB.w - 30, 16);
    }

    function frame() {
      if (phase === 'observe' || phase === 'solving') {
        const tol = phase === 'observe' ? TOL / 10 : TOL;
        const r = body.relax(0, -G, { maxSteps: 1200, tol });
        el.status.textContent = `solving update ${history.length} · residual ${r.residual.toExponential(1)}`;
        if (r.residual < tol) {
          if (phase === 'observe') {
            xObs = Float64Array.from(body.x);
            Xu = Float64Array.from(xObs);
            let peak = 0;
            for (let e = 0; e < mesh.m; e++) peak = Math.max(peak, body.mises[e]);
            vmax = peak * 1.6;
            buildLegend(el.legend, vmax, 'Mises stress');
            // Panel (b) must always show the current guess under gravity, so solve
            // the naive guess (the scan itself) right away.
            say('The red dashes in (b) are the scan: tissue <b>already deformed by gravity</b>. ' +
                'Our first guess in (a) is the naive one, the scan itself. Applying gravity to it…');
            beginIteration();
          } else {
            finishIteration();
          }
        }
      } else if (phase === 'ready' && auto) {
        if (--budget <= 0) beginIteration();
      } else {
        el.status.textContent = phase === 'done'
          ? `converged in ${history.length - 1} updates · α = ${alpha.toFixed(2)}`
          : `${mesh.n} nodes · ${mesh.m} linear triangles`;
      }
      draw();
    }
    viewA.onResize = viewB.onResize = () => { draw(); drawChart(); };
    start();
    whileVisible(canvasB, frame);
  }

  // ---------------------------------------------------------------------------
  // Figure 2: configuration-mapping diagram
  // ---------------------------------------------------------------------------

  const CONFIG_INFO = {
    ref: ['Reference domain',
      'A fixed computational mesh. It never moves during the optimization, so the standard finite element machinery stays untouched.',
      'X ∈ Ω_X'],
    unl: ['Unloaded configuration',
      'The unknown stress-free shape: what the tissue would look like without gravity or contact. It is described as a displacement field over the reference mesh and is one of the quantities we optimize.',
      'X̃ = X + u₀'],
    pred: ['Predicted loaded configuration',
      'Where the current guess ends up once the known loads, such as gravity, act on it. It comes from a full nonlinear hyperelastic finite element solve.',
      'x = X + u'],
    obs: ['Observed configurations',
      'Shapes measured from imaging, for example MRI of the same patient in different positions. Using two or more is what lets the method tell geometry and stiffness apart.',
      'x̂ from imaging'],
    u0: ['Geometry unknown, u₀',
      'Maps the reference mesh to the unloaded shape. Its deformation gradient enters the stress calculation, so gradients flow back to the geometry.',
      'F₀ = ∇u₀ + I'],
    u: ['Forward solve, u',
      'The total deformation from the reference mesh to the loaded state, solved with Newton’s method inside JAX-FEM.',
      'F = ∇u + I'],
    ut: ['Physical deformation, ũ',
      'The deformation measured from the stress-free state. Stress follows from a hyperelastic strain energy whose material parameters are optimized along with the geometry.',
      'F̃ = F F₀⁻¹,  P̃ = ∂W/∂F̃'],
    uh: ['Mismatch, û',
      'The difference between prediction and observation. The objective combines a nodal-position misfit with a deformation-gradient misfit and is minimized using end-to-end derivatives of the discretized nonlinear FE model.',
      'J = J_x + J_F'],
  };

  function initConfigDiagram() {
    const svg = document.getElementById('config-svg');
    if (!svg) return;
    const box = document.getElementById('config-explain');
    const items = svg.querySelectorAll('[data-k]');
    const show = (k) => {
      items.forEach((n) => n.classList.toggle('hot', n.dataset.k === k));
      const [title, body, eq] = CONFIG_INFO[k];
      box.innerHTML = `<h4>${title}</h4><p>${body}</p><span class="eq">${eq}</span>`;
    };
    items.forEach((n) => {
      n.addEventListener('mouseenter', () => show(n.dataset.k));
      n.addEventListener('click', () => show(n.dataset.k));
      n.addEventListener('focus', () => show(n.dataset.k));
    });
    show('unl');
  }

  // ---------------------------------------------------------------------------
  // Navigation
  // ---------------------------------------------------------------------------

  function initNav() {
    // Mark the nav entry of the section currently in view.
    const links = document.querySelectorAll('.nav-links a');
    const map = new Map([...links].map((a) => [a.getAttribute('href').slice(1), a]));
    const spy = new IntersectionObserver((entries) => {
      entries.forEach((en) => {
        if (en.isIntersecting && map.has(en.target.id)) {
          links.forEach((a) => a.removeAttribute('aria-current'));
          map.get(en.target.id).setAttribute('aria-current', 'true');
        }
      });
    }, { rootMargin: '-40% 0px -55% 0px' });
    document.querySelectorAll('section[id]').forEach((s) => spy.observe(s));
    const year = document.getElementById('year');
    if (year) year.textContent = new Date().getFullYear();
  }

  document.addEventListener('DOMContentLoaded', () => {
    initNav();
    initHero();
    initDemo();
    initConfigDiagram();
  });
})();
