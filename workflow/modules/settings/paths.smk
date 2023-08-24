import pathlib

GPFS_BASE_HILBERT = pathlib.Path("/gpfs/project/projects/medbioinf")

TOP_ROOT_DIR = GPFS_BASE_HILBERT

WORKDIR_EVAL = TOP_ROOT_DIR.joinpath("projects/assemblies/hybrids/eval/wd")
WORKDIR_ASSEMBLY = TOP_ROOT_DIR.joinpath("projects/assemblies/hybrids/verkko/wd")

# project repo
DIR_SNAKEFILE = pathlib.Path(workflow.basedir).resolve(strict=True)

DIR_SCRIPTS = DIR_SNAKEFILE.joinpath("scripts").resolve(strict=False)

DIR_ENVS = DIR_SNAKEFILE.joinpath("envs").resolve(strict=False)

DIR_PROC = pathlib.Path("proc")
DIR_RES = pathlib.Path("results")
DIR_LOG = pathlib.Path("log")
