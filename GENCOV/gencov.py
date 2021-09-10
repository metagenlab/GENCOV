import logging
import os
import sys
import click
import subprocess
from GENCOV import __version__

logging.basicConfig(
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M",
    format="[%(asctime)s %(levelname)s] %(message)s",
)


def get_snakefile():
    thisdir = os.path.abspath(os.path.dirname(__file__))
    sf = os.path.join(thisdir, 'Snakefile')
    if not os.path.exists(sf):
        sys.exit(f"Unable to locate the Snakemake workflow file at {sf}")
    return sf


@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    GENCOV is a snakemake workflow used to perform SARS-CoV-2 SNP-calling and consensus generation from Illumina paired-ends sequences 
    """


@cli.command(
    "run",
    context_settings=dict(ignore_unknown_options=True),
    short_help="run all GENCOV pipeline"
)
@click.option(
    "-c",
    "--configfile",
    type=click.Path(exists=True, resolve_path=True),
    help="config file",
)
@click.option(
    "-a",
    "--singularity-prefix",
    type=str,
    help="path to singularity envs",
)
@click.option(
    "-s",
    "--singularity-args",
    type=str,
    help="Singularity args (e.g -B /data:/data)",
)
@click.option(
    "-n",
    "--dryrun",
    is_flag=True,
    default=False,
    show_default=True,
    help="Snakemake dryrun to see the scheduling plan",
)
@click.option(
    "-c",
    "--cores",
    type=int,
    help="Max number of cores",
    default=2,
)


@click.argument("snakemake_args", nargs=-1, type=click.UNPROCESSED)

def run_workflow(configfile,
                 singularity_prefix,
                 singularity_args,
                 dryrun, 
                 cores,
                 snakemake_args):
    
    if dryrun:
        dryrun = '-n'
    else:
        dryrun = ''
    if singularity_prefix is None:
        singularity_prefix = os.environ['SINGULARITY_PREFIX']
    if not os.path.exists(singularity_prefix):
        logging.critical(f"conda env path not found: {singularity_prefix}")
        sys.exit(1)
    cmd = (
        f"snakemake --snakefile {get_snakefile()} --use-singularity --singularity-prefix {singularity_prefix} "
        f' --singularity-args "{singularity_args}" '
        f" --configfile {configfile} "
        f" --cores {cores} "
        f" {dryrun} "
        f"{' '.join(snakemake_args)}")
    logging.info("Executing: %s" % cmd)
    try:
        subprocess.check_call(cmd, shell=True)
    except subprocess.CalledProcessError as e:
        # removes the traceback
        logging.critical(e)
        exit(1)


if __name__ == "__main__":
    cli()