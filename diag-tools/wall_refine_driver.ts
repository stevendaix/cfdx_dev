// ============================================================
// Wall-refinement campaign driver (VMFL036, audit plan item 4).
// Executed by GitHub Actions ONLY (never in the Vibe sandbox):
//   cat diag-tools/port4_rhie_chow.ts diag-tools/wall_refine_driver.ts > campaign_wall.ts
//   npx -y tsx campaign_wall.ts
// This file is concatenated AFTER the twin, so buildMesh / makeSolver /
// campaignControls are already in scope. Controls UNCHANGED (maxIt 1500,
// sweeps 8/400, gates mom 2e-6 / pcorr 2e-7 / cont 2e-6, alphaU 0.65,
// alphaP 0.25, omega 1.8) - no gate or tolerance is modified here.
// Question under test: does Cd converge to the Tabata-Itakura reference
// (Re=50: 1.5785, Re=100: 1.0895) as the first radial ring shrinks below
// ~5% of R (near-wall sampling dominant), or does Cd stay near ~1.19
// (scheme bias in drag integration / gradient closure dominant)?
// Selfcheck before campaign delegation: re100 coarse 24x36 reproduced
// bit-exact (it=75, Cd=0.8716507056286522, relErr=0).
// ============================================================

const fsmod = require('fs');

const R = 0.5;
const ROUT = 50.0;
const CASES = [
  { tag: 're100_nr96', nr: 96, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 },
  { tag: 're100_nr72', nr: 72, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 },
  { tag: 're50_nr88', nr: 88, nt: 54, rho: 1.0, mu: 0.02, refCd: 1.5785 },
  { tag: 're100_nr128', nr: 128, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 }
];
const CHUNK_MS = 20 * 60 * 1000;
const CASE_CAP_MS = 45 * 60 * 1000;

function firstRingFrac(nr: number): number {
  return (ROUT - R) * Math.pow(1.0 / nr, 1.65) / R;
}

async function main(): Promise<void> {
  const out: any = {
    campaign: 'vmfl036_wall_refinement',
    date: '2026-10-16',
    twin: 'diag-tools/port4_rhie_chow.ts @ ci/long-campaign-diagnostics',
    controls: campaignControls(),
    question: 'Cd -> reference as first ring <= 5% R (near-wall sampling) vs Cd ~ 1.19 (scheme bias)',
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
