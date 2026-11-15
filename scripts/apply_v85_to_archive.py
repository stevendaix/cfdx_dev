"""Safely update a VERIFIED CFDX V6/V82/V83/V84 *source ZIP* to V85.

The original ZIP is never overwritten. Rebuild the output only if its entire
canonical source subtree agrees with the expected V85 manifest. Reject a
full delivery bundle containing multiple duplicate source trees: first extract
and upgrade its CF*_SOURCES_COMPLETES.zip member instead.

Requires Python 3.9+ and `git` on PATH. No extracted code is executed.
"""
from __future__ import annotations
import argparse,hashlib,json,os,shutil,stat,subprocess,tempfile,zipfile
from pathlib import Path,PurePosixPath

# Exact source version fingerprints are generated and frozen with this package.
FINGERPRINTS = {'V6': {'CMakeLists.txt': 'af8a357602a98ffcb9bddfc5acb68b6852b81512766c56e45233390bb4018152', 'apps/production_solver/main.cpp': '78548d5e1b52458f2ed27b0eb993819773e639e1d0623ff709b8b8bc2c6e6ceb', 'src/cfdx/physics/steady_incompressible_solver.h': '906df187212df85691173d572934cb4842de63095a2133c251a9a989b785b815'}, 'V82': {'CMakeLists.txt': 'e82c88050a00a320d0e03c201c752c20c2fd22042980cf41efa67fe905977128', 'apps/production_solver/main.cpp': 'd4920726bb6538df4086d886ad1e096764243429cebb42e5642838df52f2cbe7', 'src/cfdx/physics/steady_incompressible_solver.h': '7f3f05f14de4a64b272e319302a31408190193b7c3e34c1010b01f471a794d30'}, 'V83': {'CMakeLists.txt': 'b571e95fe1d5279de92d50b92eee5be9ca976c5687d1aa789eb11f4114ddc66a', 'apps/production_solver/main.cpp': 'cbb0f667fce5abd2b5422bc9ac01af333f89a1de05dd3f53ab10af026dcae95b', 'src/cfdx/physics/steady_incompressible_solver.h': '3587e9442557e0dc32096afd4d68fa917a6ec7b4249df39766059ac71f0310da'}, 'V84': {'CMakeLists.txt': 'a2b29fb09078d130835168c29727b16f775d76d91eccfd1653a419025799861d', 'apps/production_solver/main.cpp': '1a4b786db3ecd630525a558e522689cb643ed803627ac11113ee12736d8b33be', 'src/cfdx/physics/steady_incompressible_solver.h': '0c59b3154c1ce60b4d0a57b80da06bf1f6d4e1b85ceaef3e9ea3cbbbb4718f00'}}
FINGERPRINT_FILES=(
    'CMakeLists.txt',
    'apps/production_solver/main.cpp',
    'src/cfdx/physics/steady_incompressible_solver.h',
)
VERSIONS = {
    'V6':['CFDX_V85_PATCH_CUMULATIF_BFS_V6_V85.patch'],
    'V82':['CFDX_V85_PATCH_DIRECT_V82_ZIP_V85.patch'],
    'V83':['CFDX_V84_PATCH_INCREMENTAL_V83_V84.patch',
           'CFDX_V85_PATCH_INCREMENTAL_V84_V85.patch'],
    'V84':['CFDX_V85_PATCH_INCREMENTAL_V84_V85.patch'],
}
MAX_FILES=15000
MAX_UNCOMPRESSED=1200*1024*1024

def checksum(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def safe_entries(z):
    infos=z.infolist()
    if len(infos)>MAX_FILES or sum(i.file_size for i in infos)>MAX_UNCOMPRESSED:
        raise ValueError('ZIP size or member-count guard exceeded')
    seen=set()
    for info in infos:
        name=info.filename
        if '\\' in name or '\0' in name or ':' in name or name.startswith('/'):
            raise ValueError('Unsafe ZIP member: '+repr(name))
        p=PurePosixPath(name)
        if any(x in ('.git','..') for x in p.parts):
            raise ValueError('Unsafe ZIP traversal/Git metadata: '+repr(name))
        mode=(info.external_attr>>16)&0xFFFF
        if stat.S_ISLNK(mode):raise ValueError('Refuse ZIP symlink: '+repr(name))
        if name in seen:raise ValueError('Duplicate member: '+repr(name))
        seen.add(name)
    return infos

def candidate_prefixes(names):
    target='src/cfdx/physics/steady_incompressible_solver.h'
    prefixes=[]
    for name in names:
        if name.endswith(target):
            root=name[:-len(target)]
            if root+'CMakeLists.txt' in names and root+'apps/production_solver/main.cpp' in names:
                prefixes.append(root)
    return sorted(set(prefixes))

def validate_version(root):
    for version,digests in FINGERPRINTS.items():
        if digests and all((root/rel).is_file() and checksum(root/rel)==digests[rel]
                           for rel in FINGERPRINT_FILES):
            return version
    raise ValueError('Unknown/modified CFDX source preimage. No files changed. '
                     'This updater only accepts exact V6/V82/V83/V84 source versions.')

def run(input_zip,output_zip,bundle):
    input_zip,output_zip,bundle=map(Path,(input_zip,output_zip,bundle))
    if input_zip.resolve()==output_zip.resolve():
        raise ValueError('Refuse to overwrite the original archive')
    if not input_zip.is_file():raise FileNotFoundError(input_zip)
    if not bundle.is_dir():raise FileNotFoundError(bundle)
    if not shutil.which('git'):raise RuntimeError('git executable not on PATH')
    with tempfile.TemporaryDirectory(prefix='cfdx_v85_update_') as tmp:
        work=Path(tmp);dst=work/'unpacked';dst.mkdir()
        with zipfile.ZipFile(input_zip) as z:
            infos=safe_entries(z)
            names={i.filename for i in infos}
            prefixes=candidate_prefixes(names)
            if len(prefixes)!=1:
                raise ValueError('Expected exactly one source tree; this may be a full '
                                 'delivery ZIP. Upgrade its *_SOURCES_COMPLETES.zip instead.')
            prefix=prefixes[0]
            if any(name.lower().endswith('.zip') for name in names):
                raise ValueError('Nested archives detected. Upgrade the source-only ZIP '
                                 'to avoid stale duplicate source trees/reports.')
            for info in infos:
                path=dst.joinpath(*PurePosixPath(info.filename).parts)
                if info.is_dir():path.mkdir(parents=True,exist_ok=True);continue
                path.parent.mkdir(parents=True,exist_ok=True)
                with z.open(info) as src,path.open('wb') as out:shutil.copyfileobj(src,out)
                mode=(info.external_attr>>16)&0o777
                if mode: path.chmod(mode)
        root=dst.joinpath(*PurePosixPath(prefix).parts)
        version=validate_version(root)
        patch_names=VERSIONS[version]
        for name in patch_names:
            file=bundle/name
            if not file.is_file():raise FileNotFoundError('Missing bundled patch '+str(file))
            subprocess.run(['git','apply','--check',str(file)],cwd=root,check=True,
                           capture_output=True,text=True)
            subprocess.run(['git','apply',str(file)],cwd=root,check=True,
                           capture_output=True,text=True)
        manifest=json.loads((bundle/'CFDX_V85_SOURCE_MANIFEST.json').read_text())
        expected=manifest['files']
        for rel,info in expected.items():
            file=root/rel
            if not file.is_file() or checksum(file)!=info['sha256']:
                raise RuntimeError('Patch replay differs from V85 target at '+rel)
            # Older ZIP tools may not retain executability of unchanged files.
            # Only canonical V85 exec bits are restored, and every target hash
            # has already been certified. The original archive stays untouched.
            file.chmod((file.stat().st_mode & ~0o111) | (0o111 if info['executable'] else 0))
        extra=[str(p.relative_to(root)) for p in root.rglob('*')
               if p.is_file() and str(p.relative_to(root)) not in expected]
        if extra:raise RuntimeError('Unexpected source-tree files: '+repr(extra[:8]))
        meta={'status':'success','upgraded_from':version,'upgraded_to':'V85',
              'input_sha256':checksum(input_zip),'patches':patch_names,
              'source_files_verified':len(expected),'physics_qualified':False,
              'final_convergence_claim':False,
              'note':'A source update is not a native 100-iteration convergence result.'}
        # Repack the whole original ZIP namespace, including root prefix.
        # No in-place modification; write a temporary output and atomically replace.
        output_zip.parent.mkdir(parents=True,exist_ok=True)
        staged=output_zip.parent/(output_zip.name+'.partial')
        if staged.exists():staged.unlink()
        try:
            with zipfile.ZipFile(staged,'w',compression=zipfile.ZIP_DEFLATED,
                                 compresslevel=6,allowZip64=True) as z:
                for path in sorted(dst.rglob('*')):
                    if not path.is_file():continue
                    arc=path.relative_to(dst).as_posix()
                    z.write(path,arcname=arc)
                z.writestr('CFDX_V85_UPDATE_VERIFICATION.json',json.dumps(meta,indent=2))
            with zipfile.ZipFile(staged) as z:
                if z.testzip() is not None:raise RuntimeError('Output ZIP CRC invalid')
            os.replace(staged,output_zip)
        finally:
            if staged.exists():staged.unlink()
        meta['output']=str(output_zip)
        meta['output_sha256']=checksum(output_zip)
        return meta

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path)
    p.add_argument('--bundle',type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args()
    print(json.dumps(run(a.input,a.output,a.bundle),indent=2))