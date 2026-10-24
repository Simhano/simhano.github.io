/*
 * fem.js — a tiny 2D finite element toolkit that runs in the browser.
 *
 *   - Delaunay meshing (Bowyer–Watson) of arbitrary domains with holes
 *   - Linear triangles (CST), compressible neo-Hookean, plane strain
 *   - Explicit dynamics (with Kelvin–Voigt damping) and dynamic relaxation
 *   - Turbo colormap with Abaqus-style discrete contour bands
 */
(function (global) {
  'use strict';

  // ---------------------------------------------------------------------------
  // Meshing
  // ---------------------------------------------------------------------------

  function circumcircle(P, a, b, c) {
    const ax = P[a][0], ay = P[a][1];
    const bx = P[b][0], by = P[b][1];
    const cx = P[c][0], cy = P[c][1];
    const d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
    const a2 = ax * ax + ay * ay, b2 = bx * bx + by * by, c2 = cx * cx + cy * cy;
    const ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
    const uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
    return { a, b, c, cx: ux, cy: uy, r2: (ax - ux) ** 2 + (ay - uy) ** 2 };
  }

  /** Bowyer–Watson Delaunay triangulation. Returns [[i,j,k], ...] in CCW order. */
  function delaunay(pts) {
    const n = pts.length;
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const [x, y] of pts) {
      if (x < minX) minX = x; if (x > maxX) maxX = x;
      if (y < minY) minY = y; if (y > maxY) maxY = y;
    }
    const d = Math.max(maxX - minX, maxY - minY) * 20;
    const mx = (minX + maxX) / 2, my = (minY + maxY) / 2;
    const P = pts.concat([[mx - d, my - d], [mx + d, my - d], [mx, my + d]]);
    let tris = [circumcircle(P, n, n + 1, n + 2)];

    for (let p = 0; p < n; p++) {
      const px = P[p][0], py = P[p][1];
      const keep = [], edges = new Map();
      for (const t of tris) {
        const dx = px - t.cx, dy = py - t.cy;
        if (dx * dx + dy * dy < t.r2) {
          for (const [u, v] of [[t.a, t.b], [t.b, t.c], [t.c, t.a]]) {
            const key = u < v ? u * 1e6 + v : v * 1e6 + u;
            if (edges.has(key)) edges.delete(key); else edges.set(key, [u, v]);
          }
        } else keep.push(t);
      }
      tris = keep;
      for (const [u, v] of edges.values()) tris.push(circumcircle(P, u, v, p));
    }

    const out = [];
    for (const t of tris) {
      if (t.a >= n || t.b >= n || t.c >= n) continue;
      const [ax, ay] = P[t.a], [bx, by] = P[t.b], [cx, cy] = P[t.c];
      const area2 = (bx - ax) * (cy - ay) - (cx - ax) * (by - ay);
      if (Math.abs(area2) < 1e-12) continue;
      out.push(area2 > 0 ? [t.a, t.b, t.c] : [t.a, t.c, t.b]);
    }
    return out;
  }

  /**
   * Mesh a domain from a dense boundary point set plus an inside() test.
   * Interior nodes are placed on a jittered hexagonal lattice of spacing h.
   */
  function meshDomain({ boundary, inside, bbox, h, seed = 7 }) {
    let s = seed;
    const rand = () => ((s = (s * 16807) % 2147483647) / 2147483647 - 0.5);
    const pts = boundary.map(p => [p[0], p[1]]);
    const [x0, y0, x1, y1] = bbox;
    const dy = h * Math.sqrt(3) / 2;
    const minD2 = (0.6 * h) ** 2;
    for (let j = 0, y = y0; y <= y1; j++, y += dy) {
      for (let x = x0 + (j % 2 ? h / 2 : 0); x <= x1; x += h) {
        const px = x + rand() * 1e-3 * h, py = y + rand() * 1e-3 * h;
        if (!inside(px, py)) continue;
        let ok = true;
        for (const b of boundary) {
          if ((b[0] - px) ** 2 + (b[1] - py) ** 2 < minD2) { ok = false; break; }
        }
        if (ok) pts.push([px, py]);
      }
    }
    const tris = delaunay(pts).filter(([a, b, c]) => {
      const gx = (pts[a][0] + pts[b][0] + pts[c][0]) / 3;
      const gy = (pts[a][1] + pts[b][1] + pts[c][1]) / 3;
      return inside(gx, gy);
    });
    return buildMesh(pts, tris);
  }

  function buildMesh(pts, tris) {
    const n = pts.length, m = tris.length;
    const X = new Float64Array(2 * n);
    pts.forEach((p, i) => { X[2 * i] = p[0]; X[2 * i + 1] = p[1]; });
    const T = new Int32Array(3 * m);
    tris.forEach((t, e) => { T[3 * e] = t[0]; T[3 * e + 1] = t[1]; T[3 * e + 2] = t[2]; });

    // Unique edges and boundary edges (edges owned by a single triangle).
    const count = new Map();
    for (const [a, b, c] of tris) {
      for (const [u, v] of [[a, b], [b, c], [c, a]]) {
        const key = u < v ? u * 1e6 + v : v * 1e6 + u;
        const rec = count.get(key);
        if (rec) rec.n++; else count.set(key, { u, v, n: 1 });
      }
    }
    const edges = [], bnd = [];
    for (const r of count.values()) {
      edges.push(r.u, r.v);
      if (r.n === 1) bnd.push(r.u, r.v);
    }
    return { n, m, X, T, edges: Int32Array.from(edges), boundaryEdges: Int32Array.from(bnd) };
  }

  /** Sample a closed polyline / circle densely at spacing h. */
  function sampleCircle(cx, cy, r, h) {
    const k = Math.max(8, Math.round(2 * Math.PI * r / h));
    const out = [];
    for (let i = 0; i < k; i++) {
      const t = 2 * Math.PI * i / k;
      out.push([cx + r * Math.cos(t), cy + r * Math.sin(t)]);
    }
    return out;
  }

  function sampleSegment(ax, ay, bx, by, h, includeEnd) {
    const L = Math.hypot(bx - ax, by - ay);
    const k = Math.max(1, Math.round(L / h));
    const out = [];
    for (let i = 0; i < k + (includeEnd ? 1 : 0); i++) {
      out.push([ax + (bx - ax) * i / k, ay + (by - ay) * i / k]);
    }
    return out;
  }

  // ---------------------------------------------------------------------------
  // Neo-Hookean body (CST, plane strain)
  // ---------------------------------------------------------------------------

  class Body {
    constructor(mesh, { mu = 1, lambda = 4, rho = 1, eta = 0.02 } = {}) {
      this.mesh = mesh;
      this.mu = mu; this.lambda = lambda; this.rho = rho; this.eta = eta;
      const n = mesh.n, m = mesh.m;
      this.x = Float64Array.from(mesh.X);
      this.v = new Float64Array(2 * n);
      this.f = new Float64Array(2 * n);
      this.mass = new Float64Array(n);
      this.fixed = new Uint8Array(n);
      this.anchor = Float64Array.from(mesh.X); // prescribed positions of fixed nodes
      this.Dinv = new Float64Array(4 * m);
      this.A0 = new Float64Array(m);
      this.mises = new Float64Array(m);
      this.J = new Float64Array(m);
      this.setReference(mesh.X);
    }

    /** (Re)define the stress-free configuration that strains are measured from. */
    setReference(Xr) {
      const { T, m } = this.mesh;
      this.Xref = Float64Array.from(Xr);
      this.mass.fill(0);
      this.hmin = Infinity;
      for (let e = 0; e < m; e++) {
        const i = T[3 * e], j = T[3 * e + 1], k = T[3 * e + 2];
        const a = Xr[2 * j] - Xr[2 * i], c = Xr[2 * k] - Xr[2 * i];
        const b = Xr[2 * j + 1] - Xr[2 * i + 1], d = Xr[2 * k + 1] - Xr[2 * i + 1];
        const det = a * d - b * c; // Dm = [[a, c], [b, d]]
        this.Dinv[4 * e] = d / det; this.Dinv[4 * e + 1] = -c / det;
        this.Dinv[4 * e + 2] = -b / det; this.Dinv[4 * e + 3] = a / det;
        const A = Math.abs(det) / 2;
        this.A0[e] = A;
        const mm = this.rho * A / 3;
        this.mass[i] += mm; this.mass[j] += mm; this.mass[k] += mm;
        this.hmin = Math.min(this.hmin, 2 * A / Math.max(Math.hypot(a, b), Math.hypot(c, d)));
      }
    }

    /** Stable explicit time step for the current material. */
    stableDt(safety = 0.4) {
      const c = Math.sqrt((this.lambda + 2 * this.mu) / this.rho);
      return safety * this.hmin / c;
    }

    fix(predicate) {
      const X = this.mesh.X;
      for (let i = 0; i < this.mesh.n; i++) {
        if (predicate(X[2 * i], X[2 * i + 1], i)) this.fixed[i] = 1;
      }
    }

    /** Assemble nodal forces: internal (elastic + viscous) and gravity. */
    forces(gx, gy, withStress) {
      const { T, m, n } = this.mesh;
      const x = this.x, v = this.v, f = this.f, D = this.Dinv;
      const mu = this.mu, lam = this.lambda, eta = this.eta;
      for (let i = 0; i < n; i++) {
        f[2 * i] = this.mass[i] * gx;
        f[2 * i + 1] = this.mass[i] * gy;
      }
      for (let e = 0; e < m; e++) {
        const i = T[3 * e], j = T[3 * e + 1], k = T[3 * e + 2];
        const d00 = D[4 * e], d01 = D[4 * e + 1], d10 = D[4 * e + 2], d11 = D[4 * e + 3];
        const s00 = x[2 * j] - x[2 * i], s01 = x[2 * k] - x[2 * i];
        const s10 = x[2 * j + 1] - x[2 * i + 1], s11 = x[2 * k + 1] - x[2 * i + 1];
        const F00 = s00 * d00 + s01 * d10, F01 = s00 * d01 + s01 * d11;
        const F10 = s10 * d00 + s11 * d10, F11 = s10 * d01 + s11 * d11;
        const J = F00 * F11 - F01 * F10;
        const Jc = J > 0.05 ? J : 0.05; // guard against element inversion
        const lnJ = Math.log(Jc);
        const i00 = F11 / Jc, i01 = -F10 / Jc, i10 = -F01 / Jc, i11 = F00 / Jc; // F^-T
        let P00 = mu * (F00 - i00) + lam * lnJ * i00;
        let P01 = mu * (F01 - i01) + lam * lnJ * i01;
        let P10 = mu * (F10 - i10) + lam * lnJ * i10;
        let P11 = mu * (F11 - i11) + lam * lnJ * i11;
        if (eta > 0) {
          const w00 = v[2 * j] - v[2 * i], w01 = v[2 * k] - v[2 * i];
          const w10 = v[2 * j + 1] - v[2 * i + 1], w11 = v[2 * k + 1] - v[2 * i + 1];
          P00 += eta * (w00 * d00 + w01 * d10); P01 += eta * (w00 * d01 + w01 * d11);
          P10 += eta * (w10 * d00 + w11 * d10); P11 += eta * (w10 * d01 + w11 * d11);
        }
        const A = this.A0[e];
        const H00 = -A * (P00 * d00 + P01 * d01), H01 = -A * (P00 * d10 + P01 * d11);
        const H10 = -A * (P10 * d00 + P11 * d01), H11 = -A * (P10 * d10 + P11 * d11);
        f[2 * j] += H00; f[2 * j + 1] += H10;
        f[2 * k] += H01; f[2 * k + 1] += H11;
        f[2 * i] -= H00 + H01; f[2 * i + 1] -= H10 + H11;

        if (withStress) {
          // Cauchy stress sigma = mu/J (b - I) + lambda lnJ / J I  (plane strain)
          const b00 = F00 * F00 + F01 * F01, b11 = F10 * F10 + F11 * F11, b01 = F00 * F10 + F01 * F11;
          const p = lam * lnJ / Jc;
          const sxx = mu / Jc * (b00 - 1) + p, syy = mu / Jc * (b11 - 1) + p;
          const sxy = mu / Jc * b01, szz = p;
          this.mises[e] = Math.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3 * sxy * sxy);
          this.J[e] = J;
        }
      }
    }

    /** One semi-implicit Euler step. damping is a mass-proportional rate. */
    step(dt, { gx = 0, gy = 0, damping = 0, withStress = false, extra = null } = {}) {
      const n = this.mesh.n, x = this.x, v = this.v, f = this.f;
      this.forces(gx, gy, withStress);
      if (extra) extra(f);
      const decay = Math.exp(-damping * dt);
      let ke = 0;
      for (let i = 0; i < n; i++) {
        if (this.fixed[i]) {
          v[2 * i] = v[2 * i + 1] = 0;
          x[2 * i] = this.anchor[2 * i]; x[2 * i + 1] = this.anchor[2 * i + 1];
          continue;
        }
        const im = 1 / this.mass[i];
        v[2 * i] = (v[2 * i] + dt * f[2 * i] * im) * decay;
        v[2 * i + 1] = (v[2 * i + 1] + dt * f[2 * i + 1] * im) * decay;
        x[2 * i] += dt * v[2 * i];
        x[2 * i + 1] += dt * v[2 * i + 1];
        ke += this.mass[i] * (v[2 * i] ** 2 + v[2 * i + 1] ** 2);
      }
      return 0.5 * ke;
    }

    /**
     * Dynamic relaxation to static equilibrium using kinetic damping:
     * whenever kinetic energy peaks, velocities are zeroed.
     */
    relax(gx, gy, { maxSteps = 6000, tol = 1e-6 } = {}) {
      const dt = this.stableDt(0.5);
      const eta = this.eta; this.eta = 0;
      let prev = 0, steps = 0, fmax = Infinity;
      for (; steps < maxSteps; steps++) {
        const ke = this.step(dt, { gx, gy });
        if (ke < prev) this.v.fill(0);
        prev = ke;
        if (steps % 50 === 49) {
          fmax = this.residual(gx, gy);
          if (fmax < tol) break;
        }
      }
      this.eta = eta;
      this.forces(gx, gy, true);
      return { steps, residual: fmax };
    }

    /** Max nodal out-of-balance force per unit mass on free nodes. */
    residual(gx, gy) {
      this.forces(gx, gy, false);
      let r = 0;
      for (let i = 0; i < this.mesh.n; i++) {
        if (this.fixed[i]) continue;
        r = Math.max(r, Math.hypot(this.f[2 * i], this.f[2 * i + 1]) / this.mass[i]);
      }
      return r;
    }

    /** Area-weighted average of an element field onto the nodes. */
    nodal(field, out) {
      const { T, m, n } = this.mesh;
      const w = this._nw || (this._nw = new Float64Array(n));
      out = out || new Float64Array(n);
      out.fill(0); w.fill(0);
      for (let e = 0; e < m; e++) {
        const A = this.A0[e];
        for (let q = 0; q < 3; q++) { out[T[3 * e + q]] += A * field[e]; w[T[3 * e + q]] += A; }
      }
      for (let i = 0; i < n; i++) out[i] /= w[i] || 1;
      return out;
    }

    nearestNode(px, py, rmax) {
      let best = -1, bd = rmax * rmax;
      for (let i = 0; i < this.mesh.n; i++) {
        const d = (this.x[2 * i] - px) ** 2 + (this.x[2 * i + 1] - py) ** 2;
        if (d < bd) { bd = d; best = i; }
      }
      return best;
    }

    reset() {
      this.x.set(this.mesh.X);
      this.v.fill(0);
    }
  }

  // ---------------------------------------------------------------------------
  // Colormap (Google "Turbo" polynomial fit) with discrete contour bands
  // ---------------------------------------------------------------------------

  function turbo(t) {
    t = Math.min(1, Math.max(0, t));
    const r = 0.13572138 + t * (4.61539260 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))));
    const g = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))));
    const b = 0.10667330 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))));
    const c = v => Math.round(255 * Math.min(1, Math.max(0, v)));
    return `rgb(${c(r)},${c(g)},${c(b)})`;
  }

  /** Abaqus-like banded contour: t in [0,1] snapped to `bands` levels. */
  function banded(t, bands = 12) {
    const k = Math.min(bands - 1, Math.max(0, Math.floor(t * bands)));
    return turbo((k + 0.5) / bands);
  }

  global.FEM = { delaunay, meshDomain, buildMesh, sampleCircle, sampleSegment, Body, turbo, banded };
})(typeof window !== 'undefined' ? window : globalThis);
