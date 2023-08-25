import pathlib

GPFS_BASE_HILBERT = pathlib.Path("/gpfs/project/projects/medbioinf")

TOP_ROOT_DIR = GPFS_BASE_HILBERT

WORKDIR_EVAL = TOP_ROOT_DIR.joinpath("projects/assemblies/hybrids/eval/wd")
WORKDIR_ASSEMBLY = TOP_ROOT_DIR.joinpath("projects/assemblies/hybrids/verkko/wd")

# project repo
DIR_SNAKEFILE = pathlib.Path(workflow.basedir).resolve(strict=True)
assert DIR_SNAKEFILE.name == "workflow", DIR_SNAKEFILE

DIR_SCRIPTS = DIR_SNAKEFILE.joinpath("scripts").resolve(strict=True)

DIR_ENVS = DIR_SNAKEFILE.joinpath("envs").resolve(strict=True)

DIR_PROC = pathlib.Path("proc")
DIR_RES = pathlib.Path("results")
DIR_LOG = pathlib.Path("log")

DIR_WORKING = pathlib.Path(workflow.workdir_init).resolve(strict=True)
WORKDIR = DIR_WORKING
WD = DIR_WORKING

# prep for cluster jobs
WORKDIR.joinpath("log", "cluster_jobs", "err").mkdir(exist_ok=True, parents=True)
WORKDIR.joinpath("log", "cluster_jobs", "out").mkdir(exist_ok=True, parents=True)
