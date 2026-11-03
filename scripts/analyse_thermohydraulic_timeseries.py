"""Auditable time statistics for TH-01 / TH-07 / TH-08 probe samples.

Input is *time-resolved* CFDX probe data, not iteration history: time_s,
station,u_x_m_s,u_y_m_s,u_z_m_s[,temperature_k]. Samples must be uniformly
spaced, with every probe present at each physical time. Reports means,
unbiased Reynolds covariances, velocity/temperature RMS, turbulent heat flux,
normalized one-sided PSD (squared signal per Hz), block-mean drift and a
conservative effective-sample estimate from positive-lag autocorrelation.
The tool is an independent diagnostic, never a DNS/LES qualification oracle.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


class TimeseriesError(ValueError):
    pass


def require(ok: bool, error: str) -> None:
    if not ok:
        raise TimeseriesError(error)


def _read(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], bool]:
    series: dict[str,list[tuple[float,list[float]]]] = {}
    thermal: bool | None = None
    with path.open(newline='', encoding='utf-8') as source:
        reader = csv.DictReader(source)
        columns = {'time_s','station','u_x_m_s','u_y_m_s','u_z_m_s'}
        require(reader.fieldnames is not None and columns.issubset(reader.fieldnames),
                'missing physical sample columns')
        thermal = 'temperature_k' in reader.fieldnames
        for row in reader:
            name = row['station']
            require(bool(name) and '/' not in name and '\\' not in name,
                    'invalid station identifier')
            try:
                t = float(row['time_s'])
                values = [float(row[key]) for key in ('u_x_m_s','u_y_m_s','u_z_m_s')]
                if thermal:
                    values.append(float(row['temperature_k']))
            except (ValueError, TypeError) as exc:
                raise TimeseriesError('malformed numerical sample') from exc
            require(math.isfinite(t) and all(map(math.isfinite,values)),
                    'non-finite sample')
            if thermal:
                require(values[3]>0, 'temperature must be in absolute kelvin')
            previous = series.setdefault(name, [])
            require(not previous or t > previous[-1][0],
                    'time must increase strictly for each station')
            previous.append((t, values))
    require(bool(series), 'empty timeseries')
    times = np.asarray([v[0] for v in next(iter(series.values()))],dtype=float)
    require(times.size >= 3, 'at least three physical-time samples required')
    interval = np.diff(times)
    require(np.all(interval > 0) and
            np.max(np.abs(interval - interval[0])) <= 1e-8*interval[0],
            'nonuniform physical time step: use explicit resampling with error estimates')
    arrays: dict[str,np.ndarray] = {}
    for station, samples in sorted(series.items()):
        other = np.asarray([v[0] for v in samples])
        require(np.array_equal(times, other),
                'all stations must have exactly the same physical timestamps')
        arrays[station] = np.asarray([v[1] for v in samples], dtype=float)
    return times, arrays, bool(thermal)


def _independent_count(y: np.ndarray) -> float:
    """Positive-sequence integrated autocorrelation (no fake independent frames)."""
    n=y.size
    z=y-np.mean(y)
    variance=np.dot(z,z)
    if variance <= np.finfo(float).eps * max(1.0,np.dot(y,y)):
        return float(n)
    # Avoid O(n^2) autocorrelation, zero pad to an FFT power of two.
    size=1 << (2*n-1).bit_length()
    spectrum=np.fft.rfft(z,n=size)
    covariance=np.fft.irfft(spectrum*np.conjugate(spectrum),n=size)[:n]
    covariance /= np.arange(n,0,-1)
    correlations=covariance / covariance[0]
    summed=0.0
    for lag in range(1,min(n//2,n-1)):
        if correlations[lag] <= 0:
            break
        summed += float(correlations[lag])
    return float(np.clip(n/(1+2*summed),1,n))


def _psd(signal: np.ndarray, step: float) -> tuple[np.ndarray,np.ndarray]:
    n = len(signal)
    demeaned = signal - np.mean(signal)
    window=np.hanning(n)
    fft=np.fft.rfft(demeaned*window)
    density=(np.abs(fft)**2)/(np.sum(window**2)/step)
    if n % 2 == 0:
        density[1:-1] *= 2
    else:
        density[1:] *= 2
    frequency = np.fft.rfftfreq(n,d=step)
    return frequency,density


def compute(path: Path, *, discard: int = 0,
            min_samples: int = 64, min_effective: int = 20,
            blocks: int = 4, max_relative_drift: float = 0.15) -> dict:
    require(discard >= 0 and min_samples >= 4 and min_effective > 1 and blocks >= 2,
            'invalid sampling prerequisites')
    require(math.isfinite(max_relative_drift) and max_relative_drift > 0,
            'invalid drift tolerance')
    times, stations, thermal = _read(path)
    require(times.size-discard >= min_samples and times.size-discard >= blocks*8,
            'not enough post-spinup physical time samples')
    t = times[discard:]
    step = float(t[1]-t[0])
    n = int(t.size)
    axes=['u_x_m_s','u_y_m_s','u_z_m_s'] + (['temperature_k'] if thermal else [])
    report: dict = {
        'schema_version': 1,
        'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'time_step_s': step,
        'physical_duration_s': float(t[-1]-t[0]),
        'first_sample_time_s': float(t[0]),
        'samples': n,
        'spinup_discarded_samples': discard,
        'number_of_blocks': blocks,
        'stations': {},
        'physical_qualification': 'NOT_QUALIFIED',
        'status': 'DIAGNOSTIC_ONLY',
    }
    stationarity_all=True
    for name, values in stations.items():
        data=values[discard:]
        mean=np.mean(data, axis=0)
        centered=data-mean
        covariance=centered.T @ centered/(n-1)
        rms=np.sqrt(np.maximum(0,np.diag(covariance)))
        drift={}
        neff={}
        for j,axis in enumerate(axes):
            block_means=np.array([np.mean(section[:,j]) for section in np.array_split(data,blocks)])
            # Dimensionally consistent: normalize by mean magnitude when
            # present, otherwise fluctuation RMS at a near-zero mean station.
            scale=max(abs(float(mean[j])),float(rms[j]),1e-10)
            drift[axis]=float((np.max(block_means)-np.min(block_means))/scale)
            neff[axis]=_independent_count(data[:,j])
        stationary=all(v<=max_relative_drift for v in drift.values()) and \
            all(v>=min_effective for v in neff.values())
        stationarity_all &= stationary
        frequencies,power=_psd(data[:,0],step)
        item={
            'mean':dict(zip(axes,map(float,mean))),
            'rms':dict(zip(axes,map(float,rms))),
            'reynolds_uv_m2_s2':float(covariance[0,1]),
            'reynolds_uu_m2_s2':float(covariance[0,0]),
            'reynolds_vv_m2_s2':float(covariance[1,1]),
            'reynolds_ww_m2_s2':float(covariance[2,2]),
            'maximum_block_mean_drift':drift,
            'effective_samples_estimate':neff,
            'stationary_at_declared_gates':bool(stationary),
            'spectrum_u_x':{
                'frequencies_hz':[float(x) for x in frequencies],
                'psd_m2_s':[float(x) for x in power],
                'window':'Hann','estimator':'single_periodogram',
                'spectral_uncertainty':'NOT_ESTIMATED_SINGLE_RECORD',
            }
        }
        if thermal:
            item['turbulent_heat_flux_uT_mK_s']=float(covariance[0,3])
            item['turbulent_heat_flux_vT_mK_s']=float(covariance[1,3])
            item['temperature_variance_K2']=float(covariance[3,3])
        report['stations'][name]=item
    report['all_stations_stationary_at_declared_gates']=bool(stationarity_all)
    report['status']='STATISTICS_READY_NOT_PHYSICALLY_QUALIFIED' if stationarity_all \
        else 'INCOMPLETE_TIME_STATISTICS'
    return report


def main() -> int:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samples',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--discard-first',type=int,default=0)
    p.add_argument('--minimum-samples',type=int,default=64)
    p.add_argument('--minimum-independent-samples',type=int,default=20)
    p.add_argument('--max-block-mean-drift',type=float,default=0.15)
    p.add_argument('--require-stationary',action='store_true')
    args=p.parse_args()
    try:
        result=compute(args.samples, discard=args.discard_first,
            min_samples=args.minimum_samples,
            min_effective=args.minimum_independent_samples,
            max_relative_drift=args.max_block_mean_drift)
    except (TimeseriesError,OSError,OverflowError,ValueError) as exc:
        result={'status':'INCOMPLETE_TIME_STATISTICS',
                'physical_qualification':'NOT_QUALIFIED','error':str(exc)}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,allow_nan=False,indent=2)+'\n',encoding='utf-8')
    print(f"{result['status']} (physical qualification: NOT_QUALIFIED)")
    if 'error' in result:
        print(result['error'])
        return 2
    return 0 if not args.require_stationary or result['all_stations_stationary_at_declared_gates'] else 2


if __name__=='__main__':
    raise SystemExit(main())
