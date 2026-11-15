#!/usr/bin/env python3
"""BFS20 540: honest cell-data maps from native CFDX VTU files.

PyVista is the preferred renderer; VTK is the explicitly labeled headless
fallback. Both use VTK's original cell data (never vtkCellCenters interpolation).
These maps are diagnostics, not evidence of physical qualification.
"""
from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from PIL import Image

H_DEFAULT = 0.0127
U_IN_DEFAULT = 44.2
FIELDS = {
    "u_x": {"title": "Vitesse longitudinale", "unit": "m/s", "cmap": "turbo"},
    "u_y": {"title": "Vitesse transversale", "unit": "m/s", "cmap": "RdBu_r"},
    "u_mag": {"title": "Module de vitesse", "unit": "m/s", "cmap": "turbo"},
    "p_kpa": {"title": "Pression p exportée (relative)", "unit": "kPa", "cmap": "coolwarm"},
    "ux_deviation": {"title": "Écart Ux - vitesse d'entrée", "unit": "m/s", "cmap": "RdBu_r"},
    "k": {"title": "Énergie turbulente k", "unit": "m²/s²", "cmap": "viridis"},
    "omega_log": {"title": "log10(omega)", "unit": "log10(s⁻¹)", "cmap": "magma"},
}


def checksum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def load_vtu(path: Path, h: float = H_DEFAULT, u_in: float = U_IN_DEFAULT):
    if h <= 0 or not np.isfinite(h):
        raise ValueError("H must be finite and positive")
    if not np.isfinite(u_in):
        raise ValueError("Inlet speed must be finite")
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    raw = reader.GetOutput()
    if raw is None or raw.GetNumberOfCells() == 0:
        raise ValueError("Empty or unreadable VTU")
    required = ["u_x", "u_y", "p", "k", "omega"]
    arrays = {}
    for name in required + ["u_z", "nu_t"]:
        a = raw.GetCellData().GetArray(name)
        if a is None and name in required:
            raise ValueError("Missing required native cell field: " + name)
        if a is not None:
            v = np.asarray(vtk_to_numpy(a), dtype=np.float64)
            if v.shape != (raw.GetNumberOfCells(),) or not np.isfinite(v).all():
                raise ValueError("Non-finite/wrong shape field: " + name)
            arrays[name] = v
    ux, uy = arrays["u_x"], arrays["u_y"]
    uz = arrays.get("u_z", np.zeros_like(ux))
    arrays["u_mag"] = np.sqrt(ux**2 + uy**2 + uz**2)
    arrays["ux_deviation"] = ux - u_in
    arrays["p_kpa"] = arrays["p"] / 1000
    if np.any(arrays["omega"] <= 0):
        raise ValueError("Non-positive omega cannot be logged")
    arrays["omega_log"] = np.log10(arrays["omega"])
    for name in ["u_mag", "ux_deviation", "p_kpa", "omega_log"]:
        a = numpy_to_vtk(np.ascontiguousarray(arrays[name]), deep=True)
        a.SetName(name)
        raw.GetCellData().AddArray(a)
    pts = raw.GetPoints()
    if pts is None:
        raise ValueError("Missing mesh points")
    # Geometry is expressed as x/H, y/H, z/H in all graphics.
    transform = vtk.vtkTransform()
    transform.Scale(1.0/h, 1.0/h, 1.0/h)
    tf = vtk.vtkTransformFilter()
    tf.SetTransform(transform)
    tf.SetInputData(raw)
    tf.Update()
    scaled = vtk.vtkUnstructuredGrid()
    scaled.DeepCopy(tf.GetOutput())
    return scaled, arrays


def color_limits(arrays_list):
    # One common scale is derived from ALL requested states, with conservative
    # display anchors. Negative longitudinal velocity is never clipped.
    def lo(name): return min(float(a[name].min()) for a in arrays_list)
    def hi(name): return max(float(a[name].max()) for a in arrays_list)
    uy_span=max(13.0,abs(lo("u_y")),abs(hi("u_y")))
    deviation_span=max(10.0,abs(lo("ux_deviation")),abs(hi("ux_deviation")))
    specs = {
        "u_x": (min(-10.0,lo("u_x")), max(55.0,hi("u_x"))),
        "u_y": (-uy_span, uy_span),
        "u_mag": (0.0, max(55.0,hi("u_mag"))),
        "p_kpa": (min(-30.0,lo("p_kpa")), max(0.0,hi("p_kpa"))),
        "ux_deviation": (-deviation_span, deviation_span),
        "omega_log": (min(2.0,lo("omega_log")), max(8.0,hi("omega_log"))),
        "k": (min(0.0,lo("k")), max(0.005, hi("k"))),
    }
    # Fail closed if clipping could silently hide a physically important extreme.
    for name, (lo, hi) in specs.items():
        for arr in arrays_list:
            a = arr[name]
            if float(a.min()) < lo - 1e-7*max(1, abs(lo)) or float(a.max()) > hi + 1e-7*max(1, abs(hi)):
                raise ValueError(f"Field {name} outside fixed color scale [{lo}, {hi}], actual [{a.min()}, {a.max()}]")
    return specs


def _vtk_surface(scaled):
    import vtk
    filt = vtk.vtkDataSetSurfaceFilter()
    filt.SetInputData(scaled)
    filt.Update()
    poly = vtk.vtkPolyData()
    poly.DeepCopy(filt.GetOutput())
    if poly.GetCellData().GetArray("u_x") is None:
        raise ValueError("Surface filter lost cell arrays")
    return poly


def _frame(bounds, width, height):
    xlo, xhi, ylo, yhi = bounds
    cx, cy = (xlo+xhi)/2, (ylo+yhi)/2
    half_h = max((yhi-ylo)/2*1.05, (xhi-xlo) / (2*(width/height))*1.05)
    half_w = half_h * width/height
    return (cx-half_w, cx+half_w, cy-half_h, cy+half_h, cx, cy, half_h)


def _lut_vtk(cmap, limits):
    import vtk
    cm = matplotlib.colormaps[cmap]
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(256)
    lut.SetTableRange(*limits)
    lut.Build()
    for i in range(256):
        rgba = cm(i/255)
        lut.SetTableValue(i, *rgba)
    return lut


def render_vtk(surface, field, limits, cmap, bounds, out_path, window=(1900, 380)):
    import vtk
    width, height = window
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(surface)
    mapper.SetScalarModeToUseCellFieldData()
    mapper.SelectColorArray(field)
    mapper.SetLookupTable(_lut_vtk(cmap, limits))
    mapper.UseLookupTableScalarRangeOn()
    mapper.SetScalarRange(*limits)
    mapper.ScalarVisibilityOn()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().LightingOff()
    ren = vtk.vtkRenderer()
    ren.SetBackground(0.065, 0.095, 0.15)
    ren.AddActor(actor)
    window_vtk = vtk.vtkRenderWindow()
    window_vtk.OffScreenRenderingOn()
    window_vtk.SetSize(width, height)
    window_vtk.AddRenderer(ren)
    xmin,xmax,ymin,ymax,cx,cy,hh = _frame(bounds, width, height)
    cam=ren.GetActiveCamera()
    cam.SetFocalPoint(cx,cy,0)
    cam.SetPosition(cx,cy,100)
    cam.SetViewUp(0,1,0)
    cam.ParallelProjectionOn()
    cam.SetParallelScale(hh)
    ren.ResetCameraClippingRange()
    window_vtk.Render()
    img=vtk.vtkWindowToImageFilter()
    img.SetInput(window_vtk)
    img.Update()
    writer=vtk.vtkPNGWriter()
    writer.SetInputConnection(img.GetOutputPort())
    writer.SetFileName(str(out_path))
    writer.Write()
    window_vtk.Finalize()
    return (xmin,xmax,ymin,ymax)


def render_pyvista(scaled, field, limits, cmap, bounds, out_path, window=(1900,380)):
    import pyvista as pv
    width,height=window
    # PyVista is the requested viewer. Native cell values are unchanged.
    grid = pv.wrap(scaled)
    surface = grid.extract_surface()
    plotter=pv.Plotter(off_screen=True,window_size=window)
    plotter.set_background("#101827")
    plotter.add_mesh(surface, scalars=field, preference="cell", cmap=cmap,
                     clim=limits, show_scalar_bar=False, lighting=False)
    xmin,xmax,ymin,ymax,cx,cy,hh = _frame(bounds,width,height)
    plotter.camera_position = [(cx,cy,100),(cx,cy,0),(0,1,0)]
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = hh
    plotter.show(screenshot=str(out_path), auto_close=True)
    return (xmin,xmax,ymin,ymax)


def compose(raw, extent, field, limits, label, out, h, u_in, nneg, kind, backend):
    spec = FIELDS[field]
    fig = plt.figure(figsize=(15.5,5.6), facecolor="#0b1220")
    ax = fig.add_axes([0.06,0.24,0.88,0.53], facecolor="#101827")
    ax.imshow(Image.open(raw), extent=extent, origin="upper", interpolation="nearest")
    ax.set_xlim(extent[0],extent[1]); ax.set_ylim(extent[2],extent[3])
    ax.set_aspect("equal", adjustable="box")
    # Draw the actual step outline (not streamlines) to make geometry legible.
    ax.plot([extent[0],0,0,extent[1]], [1,1,0,0],color="#f0f6fc",lw=1.7,zorder=5)
    ax.set_xlabel("x/H — marche à x/H=0",color="#e2e8f0",fontsize=11)
    ax.set_ylabel("y/H — depuis la paroi aval",color="#e2e8f0",fontsize=11)
    ax.tick_params(colors="#cbd5e1",labelsize=9)
    for s in ax.spines.values():s.set_color("#64748b")
    fig.text(.06,.94,f"BFS 20 540  |  {spec['title']}",color="#e2e8f0",fontsize=19,weight="bold")
    fig.text(.06,.89,f"CFDX V83 · {label} · itération 1 provisoire, NON CONVERGÉE · rendu {backend.upper()}",color="#aabbd2",fontsize=10.5)
    cm=matplotlib.colormaps[spec["cmap"]]
    norm=colors.Normalize(*limits)
    cax=fig.add_axes([.19,.13,.65,.034])
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cm),cax=cax,orientation="horizontal")
    cb.ax.tick_params(colors="#e2e8f0",labelsize=10)
    cb.set_label(spec["unit"]+" · échelle fixe pour tous les états",color="#e2e8f0",fontsize=10)
    tag=f"Zone: {kind}    ·    Ux<0 : {nneg} / 20 540 cellules    ·    Uentrée={u_in:g} m/s    ·    H={h:g} m"
    fig.text(.06,.035,tag,color="#94a3b8",fontsize=9)
    fig.savefig(out,dpi=175,facecolor=fig.get_facecolor())
    plt.close(fig)


def write_sheet(images, out):
    selected=[Image.open(x).convert('RGB') for x in images]
    # 2x2 sheet, separate panels remain available in full resolution.
    target_w=1350;target_h=495
    canvas=Image.new('RGB',(target_w*2,target_h*2),(11,18,32))
    for i,im in enumerate(selected):
        im.thumbnail((target_w,target_h),Image.Resampling.LANCZOS)
        x=(i%2)*target_w+(target_w-im.width)//2
        y=(i//2)*target_h+(target_h-im.height)//2
        canvas.paste(im,(x,y))
    canvas.save(out,optimize=True)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vtu',type=Path,required=True)
    parser.add_argument('--out-dir',type=Path,required=True)
    parser.add_argument('--backend',choices=('auto','pyvista','vtk'),default='auto')
    parser.add_argument('--height',type=float,default=H_DEFAULT)
    parser.add_argument('--inlet-speed',type=float,default=U_IN_DEFAULT)
    parser.add_argument('--label',default='SIMPLE / SST provisoire')
    parser.add_argument('--other-vtu',type=Path,help='Optional second state for honest quantitative comparison')
    parser.add_argument('--fixed-scale-manifest',type=Path,
                        help='Reuse color_scales from a prior manifest; abort if a value falls outside')
    args=parser.parse_args(argv)
    args.out_dir.mkdir(parents=True,exist_ok=True)
    backend=args.backend
    if backend=='auto':
        try:
            import pyvista
        except ImportError:
            backend='vtk'
        else:
            backend='pyvista'
    if backend=='pyvista':
        try:
            import pyvista
        except ImportError as exc:
            raise SystemExit('PyVista absent: pip install pyvista ; ou utiliser --backend vtk') from exc
    mesh,arr=load_vtu(args.vtu,args.height,args.inlet_speed)
    if mesh.GetNumberOfCells()!=20540:
        raise SystemExit('Expected BFS20540; refuse to label another mesh as the qualification mesh')
    other_arr=None
    if args.other_vtu:
        _,other_arr=load_vtu(args.other_vtu,args.height,args.inlet_speed)
        if len(other_arr['u_x'])!=len(arr['u_x']):raise ValueError('Different number of cells between states')
    limits=color_limits([arr]+([other_arr] if other_arr is not None else []))
    if args.fixed_scale_manifest is not None:
        old=json.loads(args.fixed_scale_manifest.read_text())['color_scales']
        limits={k:tuple(map(float,old[k])) for k in FIELDS}
        for field,(lo,hi) in limits.items():
            for x in [arr]+([other_arr] if other_arr is not None else []):
                if float(x[field].min())<lo or float(x[field].max())>hi:
                    raise ValueError(f'Field {field} exceeds fixed scale [{lo}, {hi}]. Rerender all states with an expanded common scale.')
    surface=_vtk_surface(mesh) if backend=='vtk' else None
    zoom=(-5,25,-.2,4.0)
    overview=(-20,50,-.2,9.3)
    near_step=(-2,10,-.2,2.5)
    maps=[('u_x','01_ux_zoom',zoom),('u_y','02_uy_zoom',zoom),
          ('u_mag','03_speed_zoom',zoom),('p_kpa','04_pressure_zoom',zoom),
          ('u_mag','05_speed_overview',overview),('p_kpa','06_pressure_overview',overview),
          ('ux_deviation','07_deviation_ux_zoom',zoom),
          ('ux_deviation','08_deviation_ux_near_step',near_step),
          ('p_kpa','09_pressure_near_step',near_step),
          ('k','10_k_near_step',near_step),('omega_log','11_omega_near_step',near_step)]
    products=[]
    for field,filename,bounds in maps:
        raw=args.out_dir/('.raw_'+filename+'.png')
        func=render_vtk if backend=='vtk' else render_pyvista
        extent=func(surface if backend=='vtk' else mesh,field,limits[field],FIELDS[field]['cmap'],bounds,raw)
        png=args.out_dir/(filename+'.png')
        compose(raw,extent,field,limits[field],args.label,png,args.height,args.inlet_speed,int(np.sum(arr['u_x']<0)),
                'canal et aval' if bounds==overview else 'autour de la marche',backend)
        raw.unlink(missing_ok=True)
        products.append(png)
    sheet=args.out_dir/'00_planche_cartographies_vitesse_pression.png'
    write_sheet(products[:4],sheet)
    stats={
        'source_vtu':str(args.vtu),'sha256':checksum(args.vtu),
        'case':'BFS20540 V83','iteration':1,'status':'provisional_not_converged',
        'engine':backend,'actual_pyvista_execution':backend=='pyvista',
        'cells':int(mesh.GetNumberOfCells()),'points':int(mesh.GetNumberOfPoints()),
        'step_height_m':args.height,'inlet_speed_m_per_s':args.inlet_speed,
        'data_association':'cell_data, no vtkCellCenters',
        'negative_ux_cells':int(np.sum(arr['u_x']<0)),
        'ranges':{name:{'min':float(v.min()),'max':float(v.max()),'median':float(np.median(v))} for name,v in arr.items()},
        'color_scales':{k:list(v) for k,v in limits.items()},
        'warnings':['Provisional V83 first iteration, not a physically validated BFS solution',
                    'Pressure p is exported in Pa; not a validated Cp nor a wall pressure trace',
                    'No temporal evolution can be inferred from a single VTU',
                    'No backflow (Ux < 0) in this accepted state'],
        'maps':[x.name for x in products]+[sheet.name],
    }
    if args.other_vtu:
        stats['other_vtu']=str(args.other_vtu)
        stats['other_sha256']=checksum(args.other_vtu)
        stats['max_abs_difference_between_states']={x:float(np.max(np.abs(arr[x]-other_arr[x]))) for x in ['u_x','u_y','p','k','omega','nu_t'] if x in arr and x in other_arr}
    (args.out_dir/'cartography_manifest.json').write_text(json.dumps(stats,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'backend':backend,'images':len(products)+1,'negative_Ux':stats['negative_ux_cells'],'diff':stats.get('max_abs_difference_between_states')},indent=2))
    return 0


if __name__=='__main__':
    sys.exit(main())