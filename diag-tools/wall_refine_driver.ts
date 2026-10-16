// ============================================================
// Alpha-sensitivity campaign driver v4 (VMFL036, audit item 2:
// Rhie-Chow algebraic consistency - Majumdar relaxation treatment).
// Executed by GitHub Actions ONLY (never in the Vibe sandbox):
//   cat diag-tools/port4_rhie_chow.ts diag-tools/wall_refine_driver.ts > campaign_wall.ts
//   npx -y tsx campaign_wall.ts
// Concatenated AFTER the twin: buildMesh / makeSolver / campaignControls
// in scope. This file REPLACES the v3 wall-refinement driver on this
// branch: the v3 evidence is committed at
// diag-evidence/logs/wall_refinement_2026-10-16.json and the v3 source is
// preserved in git history (blob SHA 778d5e95).
//
// v4 question: the production RC term of #895 weights the flux with the
// RELAXED diagonals (D^RC = V/apx with apx = a_P/alpha_u), while the
// pressure-correction velocity update uses D^pc = V/a_P, i.e.
// D^RC = alpha_u * D^pc (cf. docs/validation/N9_S9_RHIE_CHOW_CONSISTENCY.md,
// Majumdar section). Property to verify empirically: the CONVERGED
// solution must be independent of alpha_u - the relaxation cancels
// between the H/A predictor and the RC pressure flux at the fixed point.
//
// Only ctrl.alphaU is varied per run. Every other control is IDENTICAL to
// campaignControls() - no gate, no tolerance, no reference is modified.
// Diagnostic only.
//
// Plumbing self-checks recorded in the JSON:
//  - re100 coarse at alphaU=0.65 must reproduce Cd=0.8716507056286522
//    bit-exact (2026-10-15 campaign value);
//  - if all alphas return identical iteration counts and trajectories,
//    the override is NOT plumbed through the twin and the campaign is
//    INCONCLUSIVE.
// ============================================================

const fsmod = require('fs');

const ALPHAS = [0.5, 0.65, 0.8, 0.95];
const CASES = [
  { tag: 're100_coarse_24x36', nr: 24, nt: 36, rho: 1.0, mu: 0.01, refCd: 1.0895,
    selfCheckCd: 0.8716507056286522 },
  { tag: 're100_fine_48x72', nr: 48, nt: 72, rho: 1.0, mu: 0.01, refCd: 1.0895 }
];
const CHUNK_MS = 20 * 60 * 1000;
const CASE_CAP_MS = 45 * 60 * 1000;

function alphaTag(alpha: number): string {
  return 'a' + String(alpha).replace('.', '_');
}

async function main(): Promise<void> {
  const out: any = {
    campaign: 'vmfl036_alpha_sensitivity',
    version: 4,
    date: '2026-10-16',
    twin: 'diag-tools/port4_rhie_chow.ts @ ci/long-campaign-diagnostics',
    question: 'v4: Majumdar relaxation consistency - is the converged Cd independent of alpha_u? (D^RC = alpha_u * D^pc, cf. N9_S9_RHIE_CHOW_CONSISTENCY.md)',
    controlsAtDefault: campaignControls(),
    controlsNote: 'identical to campaignControls() except alphaU (varied per case); gates unchanged',
    plumbingCheck: 're100_coarse at alphaU=0.65 must reproduce Cd=0.8716507056286522 bit-exact; distinct iteration counts across alphas required, else INCONCLUSIVE',
    cases: {}
  };
  for (const cs of CASES) {
    for (const alpha of ALPHAS) {
      const tag = cs.tag + '_' + alphaTag(alpha);
      const caseEntry: any = {
        alphaU: alpha, nr: cs.nr, nt: cs.nt, rho: cs.rho, mu: cs.mu, refCd: cs.refCd
      };
      out.cases[tag] = caseEntry;
      const t0 = Date.now();
      try {
        const mesh = buildMesh(cs.nr, cs.nt);
        const solver = makeSolver(mesh, cs.rho, cs.mu);
        const ctrl = campaignControls();
        ctrl.alphaU = alpha;
        let st: any = null;
        let res: any = null;
        let traj: string[] = [];
        caseEntry.cells = mesh.n;
        while (true) {
          res = solver.solveChunk(ctrl, st, CHUNK_MS, 50);
          traj = traj.concat(res.traj);
          st = res.state;
          console.log('[' + tag + '] chunk: it=' + res.it + ' converged=' + res.converged
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
        if (cs.selfCheckCd !== undefined && alpha === 0.65) {
          caseEntry.selfCheckCd = cs.selfCheckCd;
          caseEntry.selfCheckDelta = res.cd - cs.selfCheckCd;
        }
        console.log('[' + tag + '] DONE it=' + res.it + ' Cd=' + res.cd.toFixed(6)
          + ' (cdp=' + res.cd_p.toFixed(4) + ' cdv=' + res.cd_v.toFixed(4) + ')'
          + ' err=' + caseEntry.errPct.toFixed(2) + '%');
      } catch (err: any) {
        caseEntry.error = String(err);
        console.log('[' + tag + '] ERROR: ' + err);
      }
    }
  }
  // Spread summary per case (alpha-independence diagnostic, NOT a gate).
  const summary: any = {};
  for (const cs of CASES) {
    const cds: number[] = [];
    const its: number[] = [];
    for (const alpha of ALPHAS) {
      const e = out.cases[cs.tag + '_' + alphaTag(alpha)];
      if (e && e.converged && typeof e.cd === 'number') { cds.push(e.cd); its.push(e.it); }
    }
    if (cds.length >= 2) {
      const mn = Math.min.apply(null, cds);
      const mx = Math.max.apply(null, cds);
      const distinctIts = Array.from(new Set(its)).length;
      summary[cs.tag] = {
        convergedRuns: cds.length,
        cdMin: mn, cdMax: mx,
        cdSpreadPct: 100.0 * (mx - mn) / cs.refCd,
        distinctItCounts: distinctIts,
        plumbingOk: distinctIts > 1,
        alphaIndependentHeuristic: (mx - mn) / cs.refCd < 0.005 && distinctIts > 1
      };
    }
  }
  out.summary = summary;
  fsmod.mkdirSync('diag-evidence/logs', { recursive: true });
  fsmod.writeFileSync('diag-evidence/logs/alpha_sensitivity_2026-10-16.json', JSON.stringify(out, null, 2));
  console.log('WROTE diag-evidence/logs/alpha_sensitivity_2026-10-16.json');
}

main().catch((e) => { console.error('CAMPAIGN FAILED: ' + e); process.exit(1); });
