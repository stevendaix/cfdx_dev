// ============================================================
// Wall-refinement campaign driver v3 (VMFL036, audit plan items 3+4).
// Executed by GitHub Actions ONLY (never in the Vibe sandbox):
//   cat diag-tools/port4_rhie_chow.ts diag-tools/wall_refine_driver.ts > campaign_wall.ts
//   npx -y tsx campaign_wall.ts
// Concatenated AFTER the twin: buildMesh / makeSolver / campaignControls in
// scope. Controls UNCHANGED (maxIt 1500, sweeps 8/400, gates 2e-6/2e-7/2e-6,
// alphaU 0.65, alphaP 0.25, omega 1.8) - no gate or tolerance modified.
//
// v3: fixes the momentum-budget sign of v2. The correct steady identity over
// the fluid domain is   sum_boundary[ rho*u_x*(u.n) + p*nx - (tau.n)_x ]*A = 0
// (n out of the fluid domain), hence
//   D = -(convTotal + presTotal - viscTotal) / qA.
// v2 had D = convTotal - presTotal + viscTotal (sign error in the assembled
// cd_budget only; all recorded components were correct). v3 also adds the
// telescoped wall-force decomposition with the EQUATION closures:
//   wall_eq_pressure = sum_sphere p_face*nx*A  (extrapolated wall pressure,
//                      nx = stored normal pointing INTO the sphere)
//   wall_eq_visc     = sum_sphere (mu*A/bD)*ux_owner (two-point wall shear)
//   wall_eq_total    = wall_eq_pressure + wall_eq_visc  (= D by telescoping)
// compared against the drag() sampling (cell-center p, Gauss gradients).
// ============================================================

const fsmod = require('fs');

const R = 0.5;
const ROUT = 50.0;
const CASES = [
  { tag: 're100_nr96', nr: 96, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 },
  { tag: 're100_nr72', nr: 72, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 },
  { tag: 're50_nr88', nr: 88, nt: 54, rho: 1.0, mu: 0.02, refCd: 1.5785 },
  { tag: 're100_nr128', nr: 128, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 },
  { tag: 're100_nt144', nr: 48, nt: 144, rho: 1.0, mu: 0.01, refCd: 1.0895 },
  { tag: 're100_nr96_nt144', nr: 96, nt: 144, rho: 1.0, mu: 0.01, refCd: 1.0895 }
];
const CHUNK_MS = 20 * 60 * 1000;
const CASE_CAP_MS = 45 * 60 * 1000;

function firstRingFrac(nr: number): number {
  return (ROUT - R) * Math.pow(1.0 / nr, 1.65) / R;
}

// Far-field momentum-budget drag (correct sign, v3).
function budgetDrag(solver: any, rho: number, mu: number): any {
  const m = solver.mesh;
  const itn = solver.internals();
  const ux = itn.ux, ur = itn.ur, p = itn.p, phi = itn.phi;
  const nc = m.n;
  const guxX = new Float64Array(nc), guxR = new Float64Array(nc);
  const gurX = new Float64Array(nc), gurR = new Float64Array(nc);
  solver.gradients(ux, guxX, guxR, true, 1.0);
  solver.gradients(ur, gurX, gurR, true, 0.0);
  let convIn = 0.0, convOut = 0.0, presIn = 0.0, presOut = 0.0;
  let viscIn = 0.0, viscOut = 0.0, massIn = 0.0, massOut = 0.0;
  let nIn = 0, nOut = 0;
  for (let fid = 0; fid < m.nf; fid++) {
    if (m.fNb[fid] !== -1 || !m.fOuter[fid]) continue;
    const o = m.fOwner[fid];
    let nx = m.fNx[fid], nr = m.fNr[fid];
    const dotc = nx * (m.fCx[fid] - m.cx[o]) + nr * (m.fCr[fid] - m.cr[o]);
    if (dotc < 0.0) { nx = -nx; nr = -nr; }
    const A = m.fArea[fid];
    const inflow = m.fCx[fid] < 0.0;
    const uxf = inflow ? 1.0 : ux[o];
    const pf = inflow ? p[o] : 0.0;
    const mdot = phi[fid];
    const txx = 2.0 * mu * guxX[o];
    const txr = mu * (guxR[o] + gurX[o]);
    const txn_x = txx * nx + txr * nr;
    const conv = mdot * uxf;
    const pres = pf * nx * A;
    const visc = txn_x * A;
    if (inflow) { convIn += conv; presIn += pres; viscIn += visc; massIn += mdot; nIn++; }
    else { convOut += conv; presOut += pres; viscOut += visc; massOut += mdot; nOut++; }
  }
  const convTotal = convIn + convOut;
  const presTotal = presIn + presOut;
  const viscTotal = viscIn + viscOut;
  const qA = 0.5 * rho * 1.0 * (Math.PI * 0.25);
  const D = -(convTotal + presTotal - viscTotal);
  return {
    qA: qA,
    cd_budget: D / qA,
    massIn: massIn, massOut: massOut, massImbalance: massIn + massOut,
    facesIn: nIn, facesOut: nOut,
    convIn: convIn, convOut: convOut, convTotal: convTotal,
    presIn: presIn, presOut: presOut, presTotal: presTotal,
    viscIn: viscIn, viscOut: viscOut, viscTotal: viscTotal
  };
}

// Telescoped wall force with the EQUATION closures: extrapolated wall
// pressure (sphereW linear extrapolation, identical formula to the twin
// gradients() sphere branch) and two-point wall shear mu*A/bD*ux_owner.
// By the telescoped discrete momentum identity this must equal the budget D.
function wallEquationDrag(solver: any, mu: number): any {
  const m = solver.mesh;
  const itn = solver.internals();
  const ux = itn.ux, p = itn.p, bD = itn.bD;
  let pressure = 0.0, visc = 0.0, nFaces = 0, nExtrap = 0;
  for (let fid = 0; fid < m.nf; fid++) {
    if (m.fNb[fid] !== -1 || !m.fSphere[fid]) continue;
    const o = m.fOwner[fid];
    const A = m.fArea[fid];
    const nx = m.fNx[fid], nr = m.fNr[fid];
    // extrapolated wall pressure: find ring-1 same-theta neighbour of o
    let pf = p[o];
    let nb = -1;
    for (let k = 0; k < 4; k++) {
      const f2 = m.cellFaces[o * 4 + k];
      if (m.fNb[f2] === -1) continue;
      const other = (m.fOwner[f2] === o) ? m.fNb[f2] : m.fOwner[f2];
      if (m.ci[other] === 1 && m.cj[other] === m.cj[o]) { nb = other; break; }
    }
    if (nb !== -1) {
      const w = Math.hypot(m.fCx[fid] - m.cx[o], m.fCr[fid] - m.cr[o])
        / Math.max(1e-12, Math.hypot(m.cx[o] - m.cx[nb], m.cr[o] - m.cr[nb]));
      pf = p[o] + (p[o] - p[nb]) * w;
      nExtrap++;
    }
    pressure += pf * nx * A;
    visc += (mu * A / bD[fid]) * ux[o];
    nFaces++;
  }
  const total = pressure + visc;
  return {
    sphereFaces: nFaces, extrapolated: nExtrap,
    wall_eq_pressure: pressure, wall_eq_visc: visc, wall_eq_total: total
  };
}

async function main(): Promise<void> {
  const out: any = {
    campaign: 'vmfl036_wall_refinement',
    version: 3,
    date: '2026-10-16',
    twin: 'diag-tools/port4_rhie_chow.ts @ ci/long-campaign-diagnostics',
    controls: campaignControls(),
    question: 'v3: corrected momentum budget sign + telescoped wall-force decomposition (equation closures) vs drag() sampling vs reference.',
    selfcheck: 're100 coarse 24x36 pre-campaign run: converged it=75, Cd=0.8716507056286522 (bit-exact vs 2026-10-15 campaign)',
    cases: {}
  };
  for (const cs of CASES) {
    const caseEntry: any = {
      nr: cs.nr, nt: cs.nt, rho: cs.rho, mu: cs.mu, refCd: cs.refCd,
      firstRingFrac: firstRingFrac(cs.nr)
    };
    out.cases[cs.tag] = caseEntry;
    const t0 = Date.now();
    try {
      const mesh = buildMesh(cs.nr, cs.nt);
      const solver = makeSolver(mesh, cs.rho, cs.mu);
      const ctrl = campaignControls();
      let st: any = null;
      let res: any = null;
      let traj: string[] = [];
      caseEntry.cells = mesh.n;
      while (true) {
        res = solver.solveChunk(ctrl, st, CHUNK_MS, 50);
        traj = traj.concat(res.traj);
        st = res.state;
        console.log('[' + cs.tag + '] chunk: it=' + res.it + ' converged=' + res.converged
          + ' timedOut=' + res.timedOut + ' Cd=' + res.cd.toFixed(5));
        if (res.converged || res.it >= ctrl.maxIt) break;
        if (Date.now() - t0 > CASE_CAP_MS) { caseEntry.caseBudgetExhausted = true; break; }
      }
      caseEntry.it = res.it;
      caseEntry.converged = res.converged;
      caseEntry.maxedOut = !res.converged && res.it >= ctrl.maxIt;
      caseEntry.cd = res.cd;
      caseEntry.cd_p = res.cd_p;
      caseEntry.cd_v = res.cd_v;
      caseEntry.errPct = 100.0 * (res.cd - cs.refCd) / cs.refCd;
      caseEntry.mom = res.mom;
      caseEntry.cont = res.continuity;
      caseEntry.pcorr = res.pcorr;
      caseEntry.wallMs = Date.now() - t0;
      caseEntry.traj = traj;
      try {
        const b = budgetDrag(solver, cs.rho, cs.mu);
        caseEntry.budget = b;
        const qA = b.qA;
        const w = wallEquationDrag(solver, cs.mu);
        caseEntry.wallEq = w;
        caseEntry.wallEqCd = {
          cd_eq_p: -w.wall_eq_pressure / qA,
          cd_eq_v: w.wall_eq_visc / qA,
          cd_eq_total: w.wall_eq_total / qA
        };
        console.log('[' + cs.tag + '] BUDGET cd_budget=' + b.cd_budget.toFixed(5)
          + ' wallEq_total=' + (w.wall_eq_total / qA).toFixed(5)
          + ' (cdp_eq=' + (-w.wall_eq_pressure / qA).toFixed(4)
          + ' cdv_eq=' + (w.wall_eq_visc / qA).toFixed(4) + ')'
          + ' wall Cd=' + res.cd.toFixed(5)
          + ' (cdp=' + res.cd_p.toFixed(4) + ' cdv=' + res.cd_v.toFixed(4) + ')'
          + ' massImb=' + b.massImbalance.toExponential(2));
      } catch (berr: any) {
        caseEntry.budgetError = String(berr);
        console.log('[' + cs.tag + '] BUDGET ERROR: ' + berr);
      }
      console.log('[' + cs.tag + '] DONE it=' + res.it + ' Cd=' + res.cd.toFixed(5)
        + ' (cdp=' + res.cd_p.toFixed(4) + ' cdv=' + res.cd_v.toFixed(4) + ')'
        + ' err=' + caseEntry.errPct.toFixed(2) + '% firstRing='
        + (100.0 * caseEntry.firstRingFrac).toFixed(1) + '%R');
    } catch (err: any) {
      caseEntry.error = String(err);
      console.log('[' + cs.tag + '] ERROR: ' + err);
    }
  }
  fsmod.mkdirSync('diag-evidence/logs', { recursive: true });
  fsmod.writeFileSync('diag-evidence/logs/wall_refinement_2026-10-16.json', JSON.stringify(out, null, 2));
  console.log('WROTE diag-evidence/logs/wall_refinement_2026-10-16.json');
}

main().catch((e) => { console.error('CAMPAIGN FAILED: ' + e); process.exit(1); });
