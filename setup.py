import os
import glob
import sys
from sys import platform as _platform
from setuptools import setup, find_packages
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.command.build_py import build_py as _build_py
from setupext import \
    check_for_openmp, check_for_pyembree, read_embree_location, \
    get_mercurial_changeset_id, in_conda_env
from distutils.version import LooseVersion
import pkg_resources


if sys.version_info < (2, 7):
    print("yt currently requires Python version 2.7")
    print("certain features may fail unexpectedly and silently with older versions.")
    sys.exit(1)

try:
    distribute_ver = \
        LooseVersion(pkg_resources.get_distribution("distribute").version)
    if distribute_ver < LooseVersion("0.7.3"):
        print("Distribute is a legacy package obsoleted by setuptools.")
        print("We strongly recommend that you just uninstall it.")
        print("If for some reason you cannot do it, you'll need to upgrade it")
        print("to latest version before proceeding:")
        print("    pip install -U distribute")
        sys.exit(1)
except pkg_resources.DistributionNotFound:
    pass  # yay!

MAPSERVER_FILES = []
MAPSERVER_DIRS = [
    "",
    "leaflet",
    "leaflet/images"
]

for subdir in MAPSERVER_DIRS:
    dir_name = os.path.join("yt", "visualization", "mapserver", "html", subdir)
    files = []
    for ext in ["js", "html", "css", "png", "ico", "gif"]:
        files += glob.glob("%s/*.%s" % (dir_name, ext))
    MAPSERVER_FILES.append((dir_name, files))

SHADERS_DIR = os.path.join("yt", "visualization", "volume_rendering", "shaders")
SHADERS_FILES = glob.glob(os.path.join(SHADERS_DIR, "*.vertexshader")) + \
    glob.glob(os.path.join(SHADERS_DIR, "*.fragmentshader"))

VERSION = "3.3.1"

if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')


if check_for_openmp() is True:
    omp_args = ['-fopenmp']
else:
    omp_args = None

if os.name == "nt":
    std_libs = []
else:
    std_libs = ["m"]

cython_extensions = [
    Extension("yt.analysis_modules.photon_simulator.utils",
              ["yt/analysis_modules/photon_simulator/utils.pyx"],
              include_dirs=["yt/utilities/lib"]),
    Extension("yt.analysis_modules.ppv_cube.ppv_utils",
              ["yt/analysis_modules/ppv_cube/ppv_utils.pyx"],
              libraries=std_libs),
    Extension("yt.geometry.grid_visitors",
              ["yt/geometry/grid_visitors.pyx"],
              include_dirs=["yt/utilities/lib"],
              libraries=std_libs),
    Extension("yt.geometry.grid_container",
              ["yt/geometry/grid_container.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.geometry.oct_container",
              ["yt/geometry/oct_container.pyx",
               "yt/utilities/lib/tsearch.c"],
              include_dirs=["yt/utilities/lib"],
              libraries=std_libs),
    Extension("yt.geometry.oct_visitors",
              ["yt/geometry/oct_visitors.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.geometry.particle_oct_container",
              ["yt/geometry/particle_oct_container.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.geometry.selection_routines",
              ["yt/geometry/selection_routines.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.geometry.particle_deposit",
              ["yt/geometry/particle_deposit.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.geometry.particle_smooth",
              ["yt/geometry/particle_smooth.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.geometry.fake_octree",
              ["yt/geometry/fake_octree.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.utilities.spatial.ckdtree",
              ["yt/utilities/spatial/ckdtree.pyx"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs),
    Extension("yt.utilities.lib.bitarray",
              ["yt/utilities/lib/bitarray.pyx"],
              libraries=std_libs),
    Extension("yt.utilities.lib.bounding_volume_hierarchy",
              ["yt/utilities/lib/bounding_volume_hierarchy.pyx"],
              include_dirs=["yt/utilities/lib/"],
              extra_compile_args=omp_args,
              extra_link_args=omp_args,
              libraries=std_libs,
              depends=["yt/utilities/lib/mesh_triangulation.h"]),
    Extension("yt.utilities.lib.contour_finding",
              ["yt/utilities/lib/contour_finding.pyx"],
              include_dirs=["yt/utilities/lib/",
                            "yt/geometry/"],
              libraries=std_libs),
    Extension("yt.utilities.lib.geometry_utils",
              ["yt/utilities/lib/geometry_utils.pyx"],
              extra_compile_args=omp_args,
              extra_link_args=omp_args,
              libraries=std_libs),
    Extension("yt.utilities.lib.marching_cubes",
              ["yt/utilities/lib/marching_cubes.pyx",
               "yt/utilities/lib/fixed_interpolator.c"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs,
              depends=["yt/utilities/lib/fixed_interpolator.h"]),
    Extension("yt.utilities.lib.mesh_triangulation",
              ["yt/utilities/lib/mesh_triangulation.pyx"],
              depends=["yt/utilities/lib/mesh_triangulation.h"]),
    Extension("yt.utilities.lib.pixelization_routines",
              ["yt/utilities/lib/pixelization_routines.pyx",
               "yt/utilities/lib/pixelization_constants.c"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs,
              depends=["yt/utilities/lib/pixelization_constants.h"]),
    Extension("yt.utilities.lib.primitives",
              ["yt/utilities/lib/primitives.pyx"],
              libraries=std_libs),
    Extension("yt.utilities.lib.origami",
              ["yt/utilities/lib/origami.pyx",
               "yt/utilities/lib/origami_tags.c"],
              include_dirs=["yt/utilities/lib/"],
              depends=["yt/utilities/lib/origami_tags.h"]),
    Extension("yt.utilities.lib.grid_traversal",
              ["yt/utilities/lib/grid_traversal.pyx",
               "yt/utilities/lib/fixed_interpolator.c",
               "yt/utilities/lib/kdtree.c"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs,
              extra_compile_args=omp_args,
              extra_link_args=omp_args,
              depends=["yt/utilities/lib/kdtree.h",
                       "yt/utilities/lib/fixed_interpolator.h"]),
    Extension("yt.utilities.lib.image_samplers",
              ["yt/utilities/lib/image_samplers.pyx",
               "yt/utilities/lib/fixed_interpolator.c"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs,
              extra_compile_args=omp_args,
              extra_link_args=omp_args,
              depends=["yt/utilities/lib/fixed_interpolator.h"]),
    Extension("yt.utilities.lib.partitioned_grid",
              ["yt/utilities/lib/partitioned_grid.pyx",
               "yt/utilities/lib/fixed_interpolator.c"],
              include_dirs=["yt/utilities/lib/"],
              libraries=std_libs,
              depends=["yt/utilities/lib/fixed_interpolator.h"]),
    Extension("yt.utilities.lib.element_mappings",
              ["yt/utilities/lib/element_mappings.pyx"],
              libraries=std_libs),
    Extension("yt.utilities.lib.alt_ray_tracers",
              ["yt/utilities/lib/alt_ray_tracers.pyx"],
              libraries=std_libs),
]

lib_exts = [
    "particle_mesh_operations", "depth_first_octree", "fortran_reader",
    "interpolators", "misc_utilities", "basic_octree", "image_utilities",
    "points_in_volume", "quad_tree", "ray_integrators", "mesh_utilities",
    "amr_kdtools", "lenses",
]
for ext_name in lib_exts:
    cython_extensions.append(
        Extension("yt.utilities.lib.{}".format(ext_name),
                  ["yt/utilities/lib/{}.pyx".format(ext_name)],
                  libraries=std_libs))

lib_exts = ["write_array", "ragged_arrays", "line_integral_convolution"]
for ext_name in lib_exts:
    cython_extensions.append(
        Extension("yt.utilities.lib.{}".format(ext_name),
                  ["yt/utilities/lib/{}.pyx".format(ext_name)]))

extensions = [
    Extension("yt.analysis_modules.halo_finding.fof.EnzoFOF",
              ["yt/analysis_modules/halo_finding/fof/EnzoFOF.c",
               "yt/analysis_modules/halo_finding/fof/kd.c"],
              libraries=std_libs),
    Extension("yt.analysis_modules.halo_finding.hop.EnzoHop",
              glob.glob("yt/analysis_modules/halo_finding/hop/*.c")),
    Extension("yt.frontends.artio._artio_caller",
              ["yt/frontends/artio/_artio_caller.pyx"] +
              glob.glob("yt/frontends/artio/artio_headers/*.c"),
              include_dirs=["yt/frontends/artio/artio_headers/",
                            "yt/geometry/",
                            "yt/utilities/lib/"],
              depends=glob.glob("yt/frontends/artio/artio_headers/*.c")),
    Extension("yt.utilities.spatial._distance_wrap",
              glob.glob("yt/utilities/spatial/src/*.c")),
]

# EMBREE
if check_for_pyembree() is not None:
    embree_extensions = [
        Extension("yt.utilities.lib.mesh_construction",
                  ["yt/utilities/lib/mesh_construction.pyx"],
                  depends=["yt/utilities/lib/mesh_triangulation.h"]),
        Extension("yt.utilities.lib.mesh_traversal",
                  ["yt/utilities/lib/mesh_traversal.pyx"]),
        Extension("yt.utilities.lib.mesh_samplers",
                  ["yt/utilities/lib/mesh_samplers.pyx"]),
        Extension("yt.utilities.lib.mesh_intersection",
                  ["yt/utilities/lib/mesh_intersection.pyx"]),
    ]

    embree_prefix = os.path.abspath(read_embree_location())
    embree_inc_dir = [os.path.join(embree_prefix, 'include')]
    embree_lib_dir = [os.path.join(embree_prefix, 'lib')]
    if in_conda_env():
        conda_basedir = os.path.dirname(os.path.dirname(sys.executable))
        embree_inc_dir.append(os.path.join(conda_basedir, 'include'))
        embree_lib_dir.append(os.path.join(conda_basedir, 'lib'))

    if _platform == "darwin":
        embree_lib_name = "embree.2"
    else:
        embree_lib_name = "embree"

    for ext in embree_extensions:
        ext.include_dirs += embree_inc_dir
        ext.library_dirs += embree_lib_dir
        ext.language = "c++"
        ext.libraries += std_libs
        ext.libraries += [embree_lib_name]

    cython_extensions += embree_extensions

# ROCKSTAR
if os.path.exists("rockstar.cfg"):
    try:
        rd = open("rockstar.cfg").read().strip()
    except IOError:
        print("Reading Rockstar location from rockstar.cfg failed.")
        print("Please place the base directory of your")
        print("Rockstar install in rockstar.cfg and restart.")
        print("(ex: \"echo '/path/to/Rockstar-0.99' > rockstar.cfg\" )")
        sys.exit(1)

    rockstar_extdir = "yt/analysis_modules/halo_finding/rockstar"
    rockstar_extensions = [
        Extension("yt.analysis_modules.halo_finding.rockstar.rockstar_interface",
                  sources=[os.path.join(rockstar_extdir, "rockstar_interface.pyx")]),
        Extension("yt.analysis_modules.halo_finding.rockstar.rockstar_groupies",
                  sources=[os.path.join(rockstar_extdir, "rockstar_groupies.pyx")])
    ]
    for ext in rockstar_extensions:
        ext.library_dirs.append(rd)
        ext.libraries.append("rockstar")
        ext.define_macros.append(("THREADSAFE", ""))
        ext.include_dirs += [rd,
                             os.path.join(rd, "io"), os.path.join(rd, "util")]
    extensions += rockstar_extensions

if os.environ.get("GPERFTOOLS", "no").upper() != "NO":
    gpd = os.environ["GPERFTOOLS"]
    idir = os.path.join(gpd, "include")
    ldir = os.path.join(gpd, "lib")
    print(("INCLUDE AND LIB DIRS", idir, ldir))
    cython_extensions.append(
        Extension("yt.utilities.lib.perftools_wrap",
                  ["yt/utilities/lib/perftools_wrap.pyx"],
                  libraries=["profiler"],
                  library_dirs=[ldir],
                  include_dirs=[idir]))

class build_py(_build_py):
    def run(self):
        # honor the --dry-run flag
        if not self.dry_run:
            target_dir = os.path.join(self.build_lib, 'yt')
            src_dir = os.getcwd()
            changeset = get_mercurial_changeset_id(src_dir)
            self.mkpath(target_dir)
            with open(os.path.join(target_dir, '__hg_version__.py'), 'w') as fobj:
                fobj.write("hg_version = '%s'\n" % changeset)
        _build_py.run(self)

class build_ext(_build_ext):
    # subclass setuptools extension builder to avoid importing cython and numpy
    # at top level in setup.py. See http://stackoverflow.com/a/21621689/1382869
    def finalize_options(self):
        from Cython.Build import cythonize
        self.distribution.ext_modules[:] = cythonize(
                self.distribution.ext_modules)
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process
        # see http://stackoverflow.com/a/21621493/1382869
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

class sdist(_sdist):
    # subclass setuptools source distribution builder to ensure cython
    # generated C files are included in source distribution.
    # See http://stackoverflow.com/a/18418524/1382869
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize(cython_extensions)
        _sdist.run(self)

setup(
    name="yt",
    version=VERSION,
    description="An analysis and visualization toolkit for Astrophysical "
                + "simulations, focusing on Adaptive Mesh Refinement data "
                  "from Enzo, Orion, FLASH, and others.",
    classifiers=["Development Status :: 5 - Production/Stable",
                 "Environment :: Console",
                 "Intended Audience :: Science/Research",
                 "License :: OSI Approved :: BSD License",
                 "Operating System :: MacOS :: MacOS X",
                 "Operating System :: POSIX :: AIX",
                 "Operating System :: POSIX :: Linux",
                 "Programming Language :: C",
                 "Programming Language :: Python :: 2",
                 "Programming Language :: Python :: 2.7",
                 "Programming Language :: Python :: 3",
                 "Programming Language :: Python :: 3.4",
                 "Programming Language :: Python :: 3.5",
                 "Topic :: Scientific/Engineering :: Astronomy",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Visualization"],
    keywords='astronomy astrophysics visualization ' +
    'amr adaptivemeshrefinement',
    entry_points={'console_scripts': [
        'yt = yt.utilities.command_line:run_main',
    ],
        'nose.plugins.0.10': [
            'answer-testing = yt.utilities.answer_testing.framework:AnswerTesting'
    ]
    },
    packages=find_packages(),
    package_data = {'':['*.pxd']},
    setup_requires=[
        'numpy',
        'cython>=0.24',
    ],
    install_requires=[
        'matplotlib',
        'setuptools>=19.6',
        'sympy',
        'numpy',
        'IPython',
        'cython',
    ],
    cmdclass={'sdist': sdist, 'build_ext': build_ext, 'build_py': build_py},
    author="The yt project",
    author_email="yt-dev@lists.spacepope.org",
    url="http://yt-project.org/",
    license="BSD",
    zip_safe=False,
    scripts=["scripts/iyt"],
    data_files=MAPSERVER_FILES + [(SHADERS_DIR, SHADERS_FILES)],
    ext_modules=cython_extensions + extensions,
)

# This info about 'ckdtree' should be incorporated somehow...
#    setup(maintainer="SciPy Developers",
#          author="Anne Archibald",
#          maintainer_email="scipy-dev@scipy.org",
#          description="Spatial algorithms and data structures",
#          url="http://www.scipy.org",
#          license="SciPy License (BSD Style)",
#          **configuration(top_path='').todict()
#   )
