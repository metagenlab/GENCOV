#!/usr/bin/env python

from Bio import Entrez
import argparse
import time
import os
import sys
import logging
import subprocess
import pprint
Entrez.email = "r.w.kmiecinski@gmail.com"


class DownloadFailed(Exception):
    def __init__(self, msg, error_log="Something went wrong", critical=True):
        self.msg=msg
        self.critical=critical
        self.error_log=error_log
        super().__init__(self.msg)

class RsyncFailed(DownloadFailed):
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)
class WgetFailed(DownloadFailed):
    def __init__(self, *args, **kwargs):
        super().__init__(*args,**kwargs)

def main(snakemake=None, CMD=None):
    parser = constr_argparser()
    args = parser.parse_args(CMD)
    activate_logging(args.log, 
                     logging.DEBUG 
                     if args.debug 
                     else logging.INFO)
    logger.info("Starting NCBI Nucleotide Database Lookup...")
    logger.info("Commandline:\n{cmd}".format(cmd=" ".join(sys.argv if CMD is None else [__file__] +  CMD) ))
    ids = search_nuc(args)
    logger.debug("search_nuc() returned the following ids:\n{ids}".format(
            ids=pprint.pformat(ids)))
    

    ftp_prefix = get_summaries(ids)
    gff_link=f"{ftp_prefix}_genomic.gff.gz"
    fasta_link=f"{ftp_prefix}_genomic.fna.gz"
    try:
        if args.gff:
            try:
                download_file(gff_link, args.output)
            except RsyncFailed as e:
                logger.warning("Rsync had a problem: \n"
                               " {log}".format(log=e.error_log))

        try:
            download_file(fasta_link, args.output)
        except RsyncFailed as e:
            logger.warning("Rsync had a problem:"
                           " \n{msg}\ndetails:\n{log}".format(
                msg=e.msg,log=e.error_log))
            download_file(fasta_link,args.output, ftp=True)
    except DownloadFailed as e:
        logger.error("Could not download file(s):\n{msg}\n"
                     "details:\n{log}".format(msg=e.msg,log=e.error_log))
        exit(1)



def constr_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="Query for ncbi nucliotide Database")
    parser.add_argument("-o","--output",  help="Output Name of Sequence")
    parser.add_argument("--gff", action="store_true", help="Get GFF")
    parser.add_argument("--max_records", default=10)
    parser.add_argument("-l", "--log", default="ncbi_search.log")
    parser.add_argument("--debug", action="store_true" )  

    return parser

def activate_logging(log_path,  level=logging.INFO):
    global logger

    logger = logging.getLogger("search_ncbi")
    logger.setLevel(level)
    formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s: '
                '>>>\n %(message)s\n<<<',
            datefmt='%d-%m-%Y %I:%M:%S %p')
    fh=logging.FileHandler(log_path)
    fh.setLevel(level)
    ch=logging.StreamHandler()
    ch.setLevel(level)
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)




def search_nuc(args):
    id_list=[]
    count=0
    query = """("{query}"[Organism] OR {query}[All Fields]) """.format(
            query=args.query)
    logger.info("using query:\n {query}".format(query=query))
    with Entrez.esearch(db="assembly",
                        retmax=args.max_records,
                        term=query) as search_handle:
        records = Entrez.read(search_handle)
        id_list = records["IdList"]
        count = int(records["Count"])
    if count==0:
        raise RuntimeError("Query To Unspecific no Sequences Found")
    if count > 100:
        logger.warn("Query was very unspecific!!! ")
    logger.info(f"Query: {query}\nreturned {count} results\n"
           f' Only {args.max_records} will be shown in the next step'
           " and only the first one will be downloaded")
    return id_list


def get_summaries(ids):
    with Entrez.esummary(db="assembly",
                         id=",".join(ids)) as sum_handle:
        records = Entrez.read(sum_handle,validate=False)
        records = records['DocumentSummarySet']['DocumentSummary']
        logger.debug("Recieved the following records: {recs}".format(
            recs=pprint.pformat(records)))
        genbankftp =""
        refseqftp = ""
        for i,rec in enumerate(records):
            #epprint(rec)
            if i==0:
                genbankFtpDir=rec['FtpPath_GenBank']
                if genbankFtpDir:
                    genbankftp=f"{genbankFtpDir}/{os.path.basename(genbankFtpDir)}"
                refSeqFtpDir=rec['FtpPath_RefSeq']
                if refSeqFtpDir:
                    refseqftp=f"{refSeqFtpDir}/{os.path.basename(refSeqFtpDir)}"
            logger.info(f"    Record: {i}\n"
                    f'Title: {rec["Organism"]}\n'
                    f'Accession: {rec["AssemblyAccession"]}\n'
                    f'Status: {rec["AssemblyStatus"]}\n'
                    f'LastUpdate: {rec["LastUpdateDate"]}')
        if refseqftp:
            return refseqftp
        elif genbankftp:
            return genbankftp
        else:
            raise RuntimeError("Could not find Ftp Link for Query")


def download_file(_file, outname, tries=3,ftp=False):
    filename=os.path.basename(_file)
    rsync_file = _file.replace("ftp", "rsync")
    ftp_file = _file
    if outname is None:
        outname=""
    rsync_cmd = ["rsync","--progress","--copy-links","--recursive",
                 "--times", "--verbose", "--ignore-existing", 
                 rsync_file, outname]
    wget_cmd = (["wget"] + (["-O", outname ] if outname else []) + 
                    ["-N", ftp_file])
    eprint("running:", " ".join(wget_cmd if ftp else rsync_cmd))
    outlog = "Nothing to Report"
    ret = 0
    for _try in range(tries):
        with subprocess.Popen(wget_cmd if ftp else rsync_cmd,
                              stdout=subprocess.PIPE, 
                              stderr=subprocess.STDOUT) as proc:
            with subprocess.Popen(["tee", "/dev/stderr"], stdin=proc.stdout, 
                                  stdout=subprocess.PIPE) as tee_proc:
                out, err =tee_proc.communicate()
                outlog = mod_string(out) 
            ret = proc.wait()
        if ret == 0:
            break
        time.sleep(.5)

    if ret!=0:
        if not ftp:
            raise RsyncFailed(
                    "Rsync failed to run {tries} time(s):\n{cmd}".format(
                        cmd=" ".join(rsync_cmd), tries=tries),
                    error_log=outlog)
        else:
            raise WgetFailed("Wget failed to run {tries} time(s):\n{cmd}".format(
                        cmd=" ".join(wget_cmd), tries=tries),
                    error_log=outlog)



def eprint(*args, **kwargs):
    '''
    print function that prints to stderr

    :return: returns nothing
    '''
    print(*args, file=sys.stderr, **kwargs)

def epprint(*args, **kwargs):
    '''
    pprint function that prints to stderr

    :return: returns nothing
    '''
    pprint.pprint(*args, stream=sys.stderr, **kwargs)


def mod_string(string):
    if (sys.version_info > (3, 0)):
        return string.decode()
    else:
        return string

if __name__ == "__main__":
    main()
