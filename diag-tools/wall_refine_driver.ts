// ============================================================
// Wall-refinement campaign driver v2 (VMFL036, audit plan items 3+4).
// Executed by GitHub Actions ONLY (never in the Vibe sandbox):
//   cat diag-tools/port4_rhie_chow.ts diag-tools/wall_refine_driver.ts > campaign_wall.ts
//   npx -y tsx campaign_wall.ts
// Concatenated AFTER the twin: buildMesh / makeSolver / campaignControls in
// scope. Controls UNCHANGED (maxIt 1500, sweeps 8/400, gates 2e-6/2e-7/2e-6,
// alphaU 0.65, alphaP 0.25, omega 1.8) - no gate or tolerance modified.
//
// v2 additions (after v1 refuted the near-wall sampling hypothesis):
//  - Independent drag estimator: far-field momentum budget over the outer
//    boundary (Rout=50). D = sum_outer [ rho*u_x*(u.n)*A - p*nx*A + (tau.n)_x*A ]
//    with central/BC face values - it never touches the wall stresses, so it
//    discriminates "wall drag integration overestimates" (budget ~ 1.09)
//    from "the solution itself carries the bias" (budget ~ 1.26).
//  - Tangential refinement cases (nt=144) to probe scheme bias vs theta.
//
// v1 results (2026-10-16, committed by the previous run): Cd saturates at
// ~1.256 (+15%) at Re=100 as the first ring shrinks 8.5% -> 3.3% R;
// cd_p -> 0.64 (plausible), cd_v ~ 0.62 (lit ~0.45) carries the bias;
// Re=50 nr=88: 1.7487 (+10.8%). Near-wall sampling REFUTED as dominant.
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

// Independent far-field momentum-budget drag estimator. Uses only the
// converged interior fields (solver internals) and the OUTER boundary:
//   D = sum_outer [ mdot*ux_f - p_f*nx*A + (tau.n)_x*A ]
// where n points outward from the domain. Sphere and axis never appear.
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
  // Steady x-momentum balance over the whole fluid domain:
  //   sum_outer[ rho*u_x*(u.n) - p*nx + (tau.n)_x ] + (force on fluid at sphere) = 0
  // force on fluid at sphere = -D  =>  D = convTotal - presTotal + viscTotal.
  const convTotal = convIn + convOut;
  const presTotal = presIn + presOut;
  const viscTotal = viscIn + viscOut;
  const qA = 0.5 * rho * 1.0 * (Math.PI * 0.25);
  const D = convTotal - presTotal + viscTotal;
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

async function main(): Promise<void> {
  const out: any = {
    campaign: 'vmfl036_wall_refinement',
    version: 2,
    date: '2026-10-16',
    twin: 'diag-tools/port4_rhie_chow.ts @ ci/long-campaign-diagnostics',
    controls: campaignControls(),
    question: 'v1: near-wall sampling (REFUTED - Cd saturates ~1.256). v2: far-field momentum budget vs wall integration (post-processing vs solution bias) + tangential refinement nt=144.',
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
        console.log('[' + cs.tag + '] BUDGET cd_budget=' + b.cd_budget.toFixed(5)
          + ' (wall Cd=' + res.cd.toFixed(5) + ', massImbalance=' + b.massImbalance.toExponential(2) + ')');
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
