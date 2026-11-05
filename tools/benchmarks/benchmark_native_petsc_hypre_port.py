"""Paired CFDX native-AMG reference and anisotropy-port benchmark.

Neither PETSc nor HYPRE is linked or executed. Compare the old and new CFDX
methods on the SAME 5-point diffusion matrix and manufactured right-hand side.
The C++ programs report solver setup plus solve in seconds. No hard speed gate:
clock timings are not portable. Numerical failures always stop the benchmark.
"""
import argparse
import csv
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def parse_csv_lines(output, prefix, run):
    rows = []
    for line in output.splitlines():
        if not line.startswith(prefix + ','):
            continue
        fields = line.split(',')
        if prefix == 'aniso_csv':
            _, axis, ratio, n, krylov, pc, scaling, status, iters, residual, error, elapsed = fields
            if scaling != 'none' or krylov != 'fgmres' or pc not in (
                    'native_amg', 'smoothed_aggregation_amg'):
                continue
        else:
            _, nx, ny, ratio, axis, krylov, pc, status, iters, residual, error, elapsed = fields
            n = int(nx) * int(ny)
        if float(ratio) not in (1e3, 1e6):
            continue
        row = dict(run=run, size=int(n), axis=axis, ratio=float(ratio),
                   krylov=krylov, preconditioner=pc, status=int(status),
                   iterations=int(iters), residual=float(residual),
                   relative_error=float(error), elapsed_seconds=float(elapsed))
        if row['status'] != 0 or row['residual'] > 1e-9 or row['relative_error'] > 2e-3:
            raise ValueError(f'failed numerical gate: {row}')
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-binary', type=Path, required=True)
    parser.add_argument('--port-binary', type=Path, required=True)
    parser.add_argument('--nx', type=int, default=128)
    parser.add_argument('--ny', type=int, default=128)
    parser.add_argument('--runs', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--raw', type=Path)
    args = parser.parse_args()
    if not (1 <= args.runs <= 20 and 4 <= args.nx and 4 <= args.ny):
        parser.error('invalid runs/dimensions')

    rows = []
    for run in range(args.runs):
        for name, binary, prefix in [
            ('base', args.base_binary, 'aniso_csv'),
            ('port', args.port_binary, 'native_port_csv')]:
            command = [str(binary.resolve()), '--benchmark', str(args.nx), str(args.ny)]
            finished = subprocess.run(command, stdout=subprocess.PIPE,
                                      stderr=subprocess.STDOUT, text=True,
                                      timeout=90, check=True)
            entries = parse_csv_lines(finished.stdout, prefix, run)
            if not entries:
                raise RuntimeError(f'no benchmark results from {name}: {finished.stdout}')
            rows.extend(entries)
            print(f'{name} run {run + 1}/{args.runs}: {len(entries)} passing cases',
                  file=sys.stderr)

    grouped = defaultdict(list)
    for row in rows:
        grouped[(row['size'], row['axis'], row['ratio'], row['krylov'],
                 row['preconditioner'])].append(row)
    output = []
    for (size, axis, ratio, krylov, pc), points in sorted(grouped.items()):
        times = [point['elapsed_seconds'] for point in points]
        iterations = [point['iterations'] for point in points]
        output.append(dict(size=size, axis=axis, ratio=ratio,
                           krylov=krylov, preconditioner=pc,
                           repeats=len(points), iterations_median=statistics.median(iterations),
                           elapsed_seconds_median=statistics.median(times),
                           elapsed_seconds_min=min(times),
                           elapsed_seconds_max=max(times),
                           worst_physical_residual=max(p['residual'] for p in points),
                           worst_relative_error=max(p['relative_error'] for p in points)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    if args.raw:
        with args.raw.open('w', newline='') as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    print(f'{len(output)} methods/matrices with validated physical residuals')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
