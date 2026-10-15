// Port TypeScript du solveur VMFL036 axisymetrique (commit 3475826) - jumeau numerique.
// Transcrit operation par operation src/cfdx/physics/axisymmetric_vmfl036_solver.h.
// Fichier volontairement sans AUCUN caractere backslash (aller-retour transport sur).
// Maillage: radii = R + (Rout-R)*(i/nr)^1.65, theta = pi*j/nt, R=0.5, Rout=50.
// Correctifs de transcription appliques vs la premiere version du port:
//  1. gradients(): la fermeture radiale (sgr -= q*2*pi*A_mer AVANT la division) etait
//     ecrasee par l'affectation gr[c] = sgr/vol dans la version initiale. Corrige.
//  2. Stencil sphere: la boucle sautait les faces internes au lieu des faces de bord
//     (condition inversee vs C++: if(f.neighbour==(-1)) continue;). Corrige.
// Verifications: identites de fermeture du maillage exactes a 1e-13 pres
// (sum_f s*A*nx ~ 0, sum_f s*A*nr = 2*pi*A_mer). Campaign controls conformes au test:
// maxIt 1500, momSweeps 8, presSweeps 400, gates 2e-6/2e-7/2e-6, alphaU 0.65, alphaP 0.25, omega 1.8.

const PI = 3.141592653589793238462643383279502884;

function buildMesh(nr: number, nt: number): any {
  const R = 0.5, Rout = 50.0;
  const radii: number[] = [];
  for (let i = 0; i <= nr; i++) radii.push(R + (Rout - R) * Math.pow(i / nr, 1.65));
  const theta: number[] = [];
  for (let j = 0; j <= nt; j++) theta.push(PI * j / nt);
  const nxNode = new Float64Array((nr + 1) * (nt + 1));
  const nrNode = new Float64Array((nr + 1) * (nt + 1));
  const nid = (i: number, j: number) => i * (nt + 1) + j;
  for (let i = 0; i <= nr; i++) for (let j = 0; j <= nt; j++) {
    nxNode[nid(i, j)] = radii[i] * Math.cos(theta[j]);
    nrNode[nid(i, j)] = radii[i] * Math.sin(theta[j]);
  }
  const n = nr * nt;
  const cx = new Float64Array(n), cr = new Float64Array(n);
  const carea = new Float64Array(n), cvol = new Float64Array(n);
  const ci = new Int32Array(n), cj = new Int32Array(n);
  const cid = (i: number, j: number) => i * nt + j;
  const fOwner: number[] = [], fNb: number[] = [], fArea: number[] = [];
  const fNx: number[] = [], fNr: number[] = [], fCx: number[] = [], fCr: number[] = [];
  const fSphere: number[] = [], fOuter: number[] = [], fAxis: number[] = [], fDsA: number[] = [];
  const cellFacesList: number[][] = [];
  for (let c = 0; c < n; c++) cellFacesList.push([]);
  const pending = new Map<string, number>();
  for (let i = 0; i < nr; i++) for (let j = 0; j < nt; j++) {
    const c = cid(i, j);
    ci[c] = i; cj[c] = j;
    const v = [nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)];
    let cross = 0.0, mx = 0.0, mr = 0.0;
    for (let k = 0; k < 4; k++) {
      const p = k, z = (k + 1) % 4;
      const px = nxNode[v[p]], pr = nrNode[v[p]], zx = nxNode[v[z]], zr = nrNode[v[z]];
      const crv = px * zr - zx * pr;
      cross += crv; mx += (px + zx) * crv; mr += (pr + zr) * crv;
    }
    carea[c] = 0.5 * cross;
    cx[c] = mx / (6.0 * carea[c]); cr[c] = mr / (6.0 * carea[c]);
    cvol[c] = 2.0 * PI * (mr / 6.0);
    for (let e = 0; e < 4; e++) {
      const n0 = v[e], n1 = v[(e + 1) % 4];
      const key = Math.min(n0, n1) + '_' + Math.max(n0, n1);
      if (!pending.has(key)) {
        const id = fOwner.length;
        const px = nxNode[n0], pr = nrNode[n0], qx = nxNode[n1], qr = nrNode[n1];
        const fcx = 0.5 * (px + qx), fcr = 0.5 * (pr + qr);
        const dx = qx - px, dr = qr - pr;
        const ds = Math.hypot(dx, dr);
        const fnx = dr / ds, fnr = -dx / ds;
        const isAxis = Math.abs(pr) < 1e-14 && Math.abs(qr) < 1e-14;
        const area = isAxis ? 0.0 : 2.0 * PI * Math.max(0.0, fcr) * ds;
        fOwner.push(c); fNb.push(-1);
        fCx.push(fcx); fCr.push(fcr); fArea.push(area); fDsA.push(ds);
        fNx.push(fnx); fNr.push(fnr);
        fAxis.push(isAxis ? 1 : 0);
        fSphere.push((i === 0 && !isAxis) ? 1 : 0);
        fOuter.push((i + 1 === nr && !isAxis) ? 1 : 0);
        cellFacesList[c].push(id);
        pending.set(key, id);
      } else {
        const fid = pending.get(key)!;
        fNb[fid] = c;
        fSphere[fid] = 0; fOuter[fid] = 0; fAxis[fid] = 0;
        cellFacesList[c].push(fid);
        pending.delete(key);
      }
    }
  }
  for (let k = 0; k < fOwner.length; k++) {
    if (fNb[k] !== -1) {
      const o = fOwner[k], nb = fNb[k];
      const dx = cx[nb] - cx[o], dr = cr[nb] - cr[o];
      if (fNx[k] * dx + fNr[k] * dr < 0.0) { fNx[k] = -fNx[k]; fNr[k] = -fNr[k]; }
    }
  }
  let bs = 0, bo = 0, ba = 0;
  for (let k = 0; k < fOwner.length; k++) if (fNb[k] === -1) { bs += fSphere[k]; bo += fOuter[k]; ba += fAxis[k]; }
  if (bs !== nt || bo !== nt || ba !== 2 * nr) throw new Error('boundary topology mismatch: ' + bs + '/' + bo + '/' + ba);
  const cellFaces = new Int32Array(n * 4);
  for (let c = 0; c < n; c++) for (let k = 0; k < 4; k++) cellFaces[c * 4 + k] = cellFacesList[c][k];
  return {
    nr, nt, n, nf: fOwner.length, cx, cr, carea, cvol, ci, cj,
    fOwner: Int32Array.from(fOwner), fNb: Int32Array.from(fNb),
    fArea: Float64Array.from(fArea), fNx: Float64Array.from(fNx), fNr: Float64Array.from(fNr),
    fCx: Float64Array.from(fCx), fCr: Float64Array.from(fCr), fDs: Float64Array.from(fDsA),
    fSphere: Uint8Array.from(fSphere), fOuter: Uint8Array.from(fOuter), fAxis: Uint8Array.from(fAxis),
    cellFaces
  };
}

function campaignControls(): any {
  return { maxIt: 1500, momSweeps: 8, presSweeps: 400, momTol: 2e-6, presTol: 2e-7, contTol: 2e-6, alphaU: 0.65, alphaP: 0.25, omega: 1.8 };
}

function makeSolver(mesh: any, rho: number, mu: number): any {
  const m = mesh;
  const n = m.n;
  const ux = new Float64Array(n).fill(1.0);
  const ur = new Float64Array(n);
  const p = new Float64Array(n);
  const phi = new Float64Array(m.nf);
  const apx = new Float64Array(n).fill(1.0);
  const apr = new Float64Array(n).fill(1.0);
  const pc = new Float64Array(n);
  const slotF = new Int32Array(n * 4);
  const slotS = new Float64Array(n * 4);
  const slotFsgn = new Float64Array(n * 4);
  for (let c = 0; c < n; c++) for (let k = 0; k < 4; k++) {
    const fid = m.cellFaces[c * 4 + k];
    slotF[c * 4 + k] = fid;
    slotS[c * 4 + k] = (m.fOwner[fid] === c) ? 1.0 : -1.0;
    slotFsgn[c * 4 + k] = (m.fOwner[fid] === c) ? 1.0 : -1.0;
  }
  const gwSelf = new Float64Array(n * 4);
  const gwOther = new Float64Array(n * 4);
  const gOther = new Int32Array(n * 4);
  for (let c = 0; c < n; c++) for (let k = 0; k < 4; k++) {
    const fid = m.cellFaces[c * 4 + k];
    const nb = m.fNb[fid];
    if (nb !== -1) {
      const o = (m.fOwner[fid] === c) ? m.fNb[fid] : m.fOwner[fid];
      gOther[c * 4 + k] = o;
      const dc = Math.max(1e-12, Math.hypot(m.fCx[fid] - m.cx[c], m.fCr[fid] - m.cr[c]));
      const dn = Math.max(1e-12, Math.hypot(m.fCx[fid] - m.cx[o], m.fCr[fid] - m.cr[o]));
      const wc = 1.0 / dc, wn = 1.0 / dn;
      gwSelf[c * 4 + k] = wc / (wc + wn);
      gwOther[c * 4 + k] = wn / (wc + wn);
    } else gOther[c * 4 + k] = -1;
  }
  const sphereNb = new Int32Array(n).fill(-1);
  const sphereW = new Float64Array(n);
  for (let c = 0; c < n; c++) {
    if (m.ci[c] !== 0) continue;
    for (let k = 0; k < 4; k++) {
      const fid = m.cellFaces[c * 4 + k];
      // C++: if(f.neighbour==(-1)) continue; -> on saute les faces de BORD,
      // le voisin ring-1 (meme j) est trouve sur une face INTERNE.
      if (m.fNb[fid] === -1) continue;
      const o = (m.fOwner[fid] === c) ? m.fNb[fid] : m.fOwner[fid];
      if (m.ci[o] === 1 && m.cj[o] === m.cj[c]) {
        sphereNb[c] = o;
        for (let k2 = 0; k2 < 4; k2++) {
          const f2 = m.cellFaces[c * 4 + k2];
          if (m.fNb[f2] === -1 && m.fSphere[f2]) {
            sphereW[c] = Math.hypot(m.fCx[f2] - m.cx[c], m.fCr[f2] - m.cr[c])
              / Math.max(1e-12, Math.hypot(m.cx[c] - m.cx[o], m.cr[c] - m.cr[o]));
          }
        }
      }
    }
  }
  const muAd = new Float64Array(m.nf);
  const rhoAd = new Float64Array(m.nf);
  const bD = new Float64Array(m.nf);
  for (let fid = 0; fid < m.nf; fid++) {
    const o = m.fOwner[fid], nb = m.fNb[fid];
    if (nb !== -1) {
      const d = Math.hypot(m.cx[nb] - m.cx[o], m.cr[nb] - m.cr[o]);
      muAd[fid] = mu * m.fArea[fid] / d;
      rhoAd[fid] = rho * m.fArea[fid] / d;
    } else {
      const c = o;
      bD[fid] = Math.max(1e-12, Math.abs((m.fCx[fid] - m.cx[c]) * m.fNx[fid] + (m.fCr[fid] - m.cr[c]) * m.fNr[fid]));
    }
  }
  function boundaryUx(fid: number, c: number): number {
    if (m.fSphere[fid]) return 0.0;
    if (m.fOuter[fid] && m.fCx[fid] < 0.0) return 1.0;
    return ux[c];
  }
  function boundaryUr(fid: number, c: number): number {
    if (m.fSphere[fid] || m.fAxis[fid]) return 0.0;
    if (m.fOuter[fid] && m.fCx[fid] < 0.0) return 0.0;
    return ur[c];
  }
  function computeFluxes() {
    for (let fid = 0; fid < m.nf; fid++) {
      const o = m.fOwner[fid], nb = m.fNb[fid];
      let u: number, v: number;
      if (nb !== -1) { u = 0.5 * (ux[o] + ux[nb]); v = 0.5 * (ur[o] + ur[nb]); }
      else { u = boundaryUx(fid, o); v = boundaryUr(fid, o); }
      phi[fid] = rho * m.fArea[fid] * (u * m.fNx[fid] + v * m.fNr[fid]);
    }
  }
  function gradients(q: Float64Array, gx: Float64Array, gr: Float64Array, velocity: boolean, outerValue: number) {
    gx.fill(0.0); gr.fill(0.0);
    for (let c = 0; c < n; c++) {
      let sgx = 0.0, sgr = 0.0;
      for (let k = 0; k < 4; k++) {
        const kk = c * 4 + k;
        const fid = slotF[kk];
        let qf = q[c];
        if (gOther[kk] !== -1) {
          qf = q[c] * gwSelf[kk] + q[gOther[kk]] * gwOther[kk];
        } else if (velocity) {
          if (m.fSphere[fid]) qf = 0.0;
          else if (m.fOuter[fid] && m.fCx[fid] < 0.0) qf = outerValue;
          else if (m.fAxis[fid]) qf = (outerValue === 0.0 ? 0.0 : q[c]);
        } else if (m.fOuter[fid] && m.fCx[fid] >= 0.0) {
          qf = 0.0;
        } else if (m.fSphere[fid] && sphereNb[c] !== -1) {
          qf = q[c] + (q[c] - q[sphereNb[c]]) * sphereW[c];
        }
        const s = slotS[kk];
        const A = m.fArea[fid];
        sgx += s * qf * m.fNx[fid] * A;
        sgr += s * qf * m.fNr[fid] * A;
      }
      // Fermeture radiale: AVANT la division (comme en C++).
      sgr -= q[c] * 2.0 * PI * m.carea[c];
      gx[c] = sgx / m.cvol[c]; gr[c] = sgr / m.cvol[c];
    }
  }
  const nbId = new Int32Array(n * 4);
  const nbAn = new Float64Array(n * 4);
  const nbCnt = new Int32Array(n);
  const rhsX = new Float64Array(n), rhsR = new Float64Array(n);
  const gpX = new Float64Array(n), gpR = new Float64Array(n);
  function assembleMomentum(c: any) {
    nbCnt.fill(0);
    for (let cell = 0; cell < n; cell++) {
      let ax = 0.0, ar = 0.0;
      const V = m.cvol[cell];
      let bx = -gpX[cell] * V, br = -gpR[cell] * V;
      for (let k = 0; k < 4; k++) {
        const kk = cell * 4 + k;
        const fid = slotF[kk];
        const internal = m.fNb[fid] !== -1;
        const F = slotFsgn[kk] * phi[fid];
        if (internal) {
          const o = gOther[kk];
          const D = muAd[fid];
          const an = D + Math.max(-F, 0.0);
          ax += D + Math.max(F, 0.0);
          ar += D + Math.max(F, 0.0);
          const cnb = cell * 4 + nbCnt[cell];
          nbId[cnb] = o; nbAn[cnb] = an; nbCnt[cell]++;
        } else if (m.fAxis[fid]) {
          const D = mu * m.fArea[fid] / bD[fid];
          ax += F;
          ar += D + Math.max(F, 0.0);
        } else if (m.fOuter[fid] && m.fCx[fid] >= 0.0) {
          ax += F; ar += F;
        } else {
          const ubx = boundaryUx(fid, cell), ubr = boundaryUr(fid, cell);
          const D = mu * m.fArea[fid] / bD[fid];
          const ab = D + Math.max(-F, 0.0);
          ax += D + Math.max(F, 0.0);
          ar += D + Math.max(F, 0.0);
          bx += ab * ubx; br += ab * ubr;
        }
      }
      const hoop = mu * V / (m.cr[cell] * m.cr[cell]);
      ar += hoop;
      if (!(ax > 0.0) || !(ar > 0.0)) throw new Error('invalid momentum coefficients at cell ' + cell);
      const oldx = ux[cell], oldr = ur[cell];
      apx[cell] = ax / c.alphaU;
      apr[cell] = ar / c.alphaU;
      bx += (1.0 - c.alphaU) / c.alphaU * ax * oldx;
      br += (1.0 - c.alphaU) / c.alphaU * ar * oldr;
      rhsX[cell] = bx; rhsR[cell] = br;
    }
    for (let s = 0; s < c.momSweeps; s++) {
      for (let cell = 0; cell < n; cell++) {
        let sumx = rhsX[cell], sumr = rhsR[cell];
        for (let k = 0; k < nbCnt[cell]; k++) {
          const kk = cell * 4 + k;
          const nb = nbId[kk], an = nbAn[kk];
          sumx += an * ux[nb]; sumr += an * ur[nb];
        }
        const sx = sumx / apx[cell], sr = sumr / apr[cell];
        if (!isFinite(sx) || !isFinite(sr)) throw new Error('non-finite momentum update at cell ' + cell + ' sweep ' + s);
        ux[cell] = sx; ur[cell] = sr;
      }
    }
  }
  function momentumResidual(): number {
    let r = 0.0;
    for (let cell = 0; cell < n; cell++) {
      let rx = gpX[cell] * m.cvol[cell], rr = gpR[cell] * m.cvol[cell];
      for (let k = 0; k < 4; k++) {
        const kk = cell * 4 + k;
        const fid = slotF[kk];
        const F = slotFsgn[kk] * phi[fid];
        const inb = m.fNb[fid] !== -1;
        const nb = inb ? gOther[kk] : cell;
        const D = inb ? muAd[fid] : mu * m.fArea[fid] / bD[fid];
        const an = D + Math.max(-F, 0.0);
        const ubx = inb ? ux[nb] : boundaryUx(fid, cell);
        const ubr = inb ? ur[nb] : boundaryUr(fid, cell);
        rx += (D + Math.max(F, 0.0)) * ux[cell] - an * ubx;
        rr += (D + Math.max(F, 0.0)) * ur[cell] - an * ubr;
      }
      rr += mu * m.cvol[cell] * ur[cell] / Math.max(m.cr[cell] * m.cr[cell], 1e-20);
      r = Math.max(r, Math.max(Math.abs(rx), Math.abs(rr)) / (1.0 + Math.abs(ux[cell]) + Math.abs(ur[cell])));
    }
    return r;
  }
  const pNbId = new Int32Array(n * 4);
  const pNbA = new Float64Array(n * 4);
  const pCnt = new Int32Array(n);
  const pDiag = new Float64Array(n), pRhs = new Float64Array(n);
  function faceResponse(cell: number, fnx: number, fnr: number, c: any): number {
    const v = m.cvol[cell];
    return (fnx * fnx * v) / (c.alphaU * apx[cell]) + (fnr * fnr * v) / (c.alphaU * apr[cell]);
  }
  const cellScale = new Float64Array(n);
  function refreshCellScale() {
    for (let cell = 0; cell < n; cell++) {
      let s = 0.0;
      for (let k = 0; k < 4; k++) s += Math.abs(phi[slotF[cell * 4 + k]]);
      cellScale[cell] = s;
    }
  }
  function pressureCorrection(c: any, sweepsLog?: any[]): number {
    pc.fill(0.0);
    pDiag.fill(0.0);
    pCnt.fill(0);
    for (let cell = 0; cell < n; cell++) {
      let imbalance = 0.0;
      for (let k = 0; k < 4; k++) {
        const kk = cell * 4 + k;
        imbalance += slotFsgn[kk] * phi[slotF[kk]];
      }
      pRhs[cell] = -imbalance;
      for (let k = 0; k < 4; k++) {
        const kk = cell * 4 + k;
        const fid = slotF[kk];
        if (m.fNb[fid] !== -1) {
          const nb = gOther[kk];
          const df = 0.5 * (faceResponse(cell, m.fNx[fid], m.fNr[fid], c) + faceResponse(nb, m.fNx[fid], m.fNr[fid], c));
          const a = rhoAd[fid] * df;
          pDiag[cell] += a;
          const cp = cell * 4 + pCnt[cell];
          pNbId[cp] = nb; pNbA[cp] = -a; pCnt[cell]++;
        } else if (m.fOuter[fid] && m.fCx[fid] >= 0.0) {
          const a = rho * m.fArea[fid] * faceResponse(cell, m.fNx[fid], m.fNr[fid], c) / bD[fid];
          pDiag[cell] += a;
        }
      }
      if (!(pDiag[cell] > 0.0)) throw new Error('invalid pressure diagonal at cell ' + cell);
    }
    let maxcorr = Infinity;
    for (let s = 0; s < c.presSweeps; s++) {
      maxcorr = 0.0;
      for (let cell = 0; cell < n; cell++) {
        let sum = pRhs[cell];
        for (let k = 0; k < pCnt[cell]; k++) {
          const kk = cell * 4 + k;
          sum -= pNbA[kk] * pc[pNbId[kk]];
        }
        const old = pc[cell];
        const v = old + c.omega * (sum / pDiag[cell] - old);
        const d = Math.abs(v - old);
        if (d > maxcorr) maxcorr = d;
        pc[cell] = v;
      }
      if (sweepsLog) {
        let lr = 0.0;
        for (let cell = 0; cell < n; cell++) {
          let r = pDiag[cell] * pc[cell];
          for (let k = 0; k < pCnt[cell]; k++) {
            const kk = cell * 4 + k;
            r += pNbA[kk] * pc[pNbId[kk]];
          }
          r -= pRhs[cell];
          const sabs = Math.abs(r) / (1.0 + cellScale[cell]);
          if (sabs > lr) lr = sabs;
        }
        sweepsLog.push({ sweep: s + 1, maxcorr: maxcorr, resRel: lr });
      }
      if (maxcorr < c.presTol) {
        let resRel = 0.0;
        for (let cell = 0; cell < n; cell++) {
          let r = pDiag[cell] * pc[cell];
          for (let k = 0; k < pCnt[cell]; k++) {
            const kk = cell * 4 + k;
            r += pNbA[kk] * pc[pNbId[kk]];
          }
          r -= pRhs[cell];
          const sabs = Math.abs(r) / (1.0 + cellScale[cell]);
          if (sabs > resRel) resRel = sabs;
        }
        if (resRel < c.presTol) break;
      }
    }
    if (maxcorr < c.presTol) {
      let res = 0.0, diagMax = 0.0;
      for (let cell = 0; cell < n; cell++) {
        let r = pDiag[cell] * pc[cell];
        for (let k = 0; k < pCnt[cell]; k++) {
          const kk = cell * 4 + k;
          r += pNbA[kk] * pc[pNbId[kk]];
        }
        r -= pRhs[cell];
        if (Math.abs(r) > res) res = Math.abs(r);
        if (pDiag[cell] > diagMax) diagMax = pDiag[cell];
      }
      if (res > 100.0 * diagMax * maxcorr) throw new Error('pressure-correction operator or sign regression');
    }
    for (let cell = 0; cell < n; cell++) {
      let gx = 0.0, gr = 0.0;
      for (let k = 0; k < 4; k++) {
        const kk = cell * 4 + k;
        const fid = slotF[kk];
        let qf = pc[cell];
        if (gOther[kk] !== -1) {
          qf = pc[cell] * gwSelf[kk] + pc[gOther[kk]] * gwOther[kk];
        } else if (m.fOuter[fid] && m.fCx[fid] >= 0.0) {
          qf = 0.0;
        } else if (m.fSphere[fid] && sphereNb[cell] !== -1) {
          qf = pc[cell] + (pc[cell] - pc[sphereNb[cell]]) * sphereW[cell];
        }
        const s = slotS[kk];
        gx += s * qf * m.fNx[fid] * m.fArea[fid];
        gr += s * qf * m.fNr[fid] * m.fArea[fid];
      }
      gr -= pc[cell] * 2.0 * PI * m.carea[cell];
      gx /= m.cvol[cell]; gr /= m.cvol[cell];
      ux[cell] -= m.cvol[cell] / (c.alphaU * apx[cell]) * gx;
      ur[cell] -= m.cvol[cell] / (c.alphaU * apr[cell]) * gr;
      p[cell] += c.alphaP * pc[cell];
    }
    for (let fid = 0; fid < m.nf; fid++) {
      const o = m.fOwner[fid], nb = m.fNb[fid];
      if (nb !== -1) {
        const df = 0.5 * (faceResponse(o, m.fNx[fid], m.fNr[fid], c) + faceResponse(nb, m.fNx[fid], m.fNr[fid], c));
        phi[fid] += rhoAd[fid] * df * (pc[o] - pc[nb]);
      } else if (m.fOuter[fid] && m.fCx[fid] >= 0.0) {
        phi[fid] += rho * m.fArea[fid] * faceResponse(o, m.fNx[fid], m.fNr[fid], c) / bD[fid] * pc[o];
      }
    }
    return maxcorr;
  }
  function maxContinuity(): number {
    let r = 0.0;
    for (let cell = 0; cell < n; cell++) {
      let b = 0.0, scale = 0.0;
      for (let k = 0; k < 4; k++) {
        const kk = cell * 4 + k;
        b += slotFsgn[kk] * phi[slotF[kk]];
        scale += Math.abs(phi[slotF[kk]]);
      }
      const v = Math.abs(b) / (1.0 + scale);
      if (v > r) r = v;
    }
    return r;
  }
  function drag(): number[] {
    const guxX = new Float64Array(n), guxR = new Float64Array(n);
    const gurX = new Float64Array(n), gurR = new Float64Array(n);
    gradients(ux, guxX, guxR, true, 1.0);
    gradients(ur, gurX, gurR, true, 0.0);
    let fp = 0.0, fv = 0.0;
    for (let c = 0; c < n; c++) {
      if (m.ci[c] !== 0) continue;
      for (let k = 0; k < 4; k++) {
        const fid = m.cellFaces[c * 4 + k];
        if (m.fNb[fid] !== -1 || !m.fSphere[fid]) continue;
        const nx = -m.fNx[fid], nr = -m.fNr[fid];
        const tauxx = 2.0 * mu * guxX[c];
        const tauxr = mu * (guxR[c] + gurX[c]);
        const dS = 2.0 * PI * m.fCr[fid] * m.fDs[fid];
        fp += -(p[c]) * nx * dS;
        fv += (tauxx * nx + tauxr * nr) * dS;
      }
    }
    const qA = 0.5 * rho * 1.0 * (PI * 0.25);
    return [fp / qA, fv / qA, (fp + fv) / qA];
  }
  function exportState(): any {
    return { it: 0, ux: Array.from(ux), ur: Array.from(ur), p: Array.from(p), phi: Array.from(phi) };
  }
  function importState(st: any) {
    if (!st) return;
    for (let c = 0; c < n; c++) { ux[c] = st.ux[c]; ur[c] = st.ur[c]; p[c] = st.p[c]; }
    for (let f = 0; f < m.nf; f++) phi[f] = st.phi[f];
  }
  // Chunked campaign runner: reprend depuis un checkpoint, tourne jusqu'a
  // convergence, budget maxIt atteint, ou budget temps ecoule (ms).
  function solveChunk(c: any, st: any, budgetMs: number, logEvery: number): any {
    importState(st);
    const startIt = st ? (st.it || 0) : 0;
    if (!st) computeFluxes();
    const t0 = Date.now();
    let continuity = Infinity, mom = Infinity, pcorr = Infinity;
    let lastSweepsLog: any[] = [];
    let it = startIt;
    let converged = false;
    let timedOut = false;
    const traj: string[] = [];
    let lastCompleted = startIt;
    for (it = startIt + 1; it <= c.maxIt; it++) {
      gradients(p, gpX, gpR, false, 0.0);
      assembleMomentum(c);
      mom = momentumResidual();
      computeFluxes();
      refreshCellScale();
      const sweepsLog: any[] = [];
      pcorr = pressureCorrection(c, c.logSOR ? sweepsLog : undefined);
      lastSweepsLog = sweepsLog;
      continuity = maxContinuity();
      if (it % logEvery === 0 || it === startIt + 1 || it === 1) {
        const d = drag();
        const line = 'ITER=' + it + ' cont=' + continuity.toExponential(3) + ' mom=' + mom.toExponential(3) + ' pcorr=' + pcorr.toExponential(3) + ' Cd=' + d[2].toFixed(4) + ' cdp=' + d[0].toFixed(3) + ' cdv=' + d[1].toFixed(3);
        traj.push(line);
        console.log(line);
      }
      if (continuity < c.contTol && mom < c.momTol && pcorr < c.presTol) {
        converged = true;
        break;
      }
      lastCompleted = it;
      if (Date.now() - t0 > budgetMs) { timedOut = true; break; }
    }
    const state = exportState();
    state.it = lastCompleted;
    const d = drag();
    return { state, it: state.it, continuity, mom, pcorr, converged, timedOut, done: converged || state.it >= c.maxIt, cd: d[2], cd_p: d[0], cd_v: d[1], traj, sweepsLog: lastSweepsLog };
  }
  return { solveChunk, mesh: m, n, drag, gradients, pressureCorrection, momentumResidual, maxContinuity, computeFluxes, refreshCellScale, internals: function() { return { ux: ux, ur: ur, p: p, phi: phi, apx: apx, apr: apr, gpX: gpX, gpR: gpR, muAd: muAd, rhoAd: rhoAd, bD: bD }; } };
}
