# -*- coding: utf-8 -*-

from subprocess import call, run
from urllib import request

import click
from requests import get


def prs_prompt(*, download, plink):
    """
    Interactive function to download and calculate prs.

    Parameters
    ----------
    download : bool
        if true, PGS file will be downloaded given an ID
    plink : str
        the path to plink. if not default in environment.

    Returns
    -------
        a statement
    """
    if download:
        pgs_id = click.prompt('Please input PGS ID', type=str)
        pgs_file = download_pgs(prs_id=pgs_id)
        p = run(["zcat", pgs_file, "|", 'head', '-n', "15"], shell=True)
        if p.returncode != 0:
            click.echo('The PGS file could not be viewed using cmd, please view it manually.')
        click.echo('Please check the PGS file viewed and provide us with the needed columns.')
    else:
        pgs_file = click.prompt('Please provide path to PGS file', type=str)
    id_col = click.prompt('Please provide the ID column number')
    allele = click.prompt('Please provide the effect allele column number')
    weight = click.prompt('Please provide the effect weight column number')
    cols = ' '.join([str(id_col), str(allele), str(weight)])
    file_type = click.prompt('Do you have a VCF file or binary files?', type=click.Choice(['vcf', 'bfile']))
    input_file = click.prompt('Please provide the path to input file', type=str)
    confirm = click.confirm('Please be aware that variant ID in both input file and pgs file need to match.'
                            'Do you want to continue?')
    if not confirm:
        return 'Ok. You still have the PGS file (if downloaded) but the scores were not calculated.'
    output_file = click.prompt('Please provide an output file path', type=str)
    if file_type == 'vcf':
        p = call(
            plink + " --vcf " + input_file + " --score " + pgs_file + " " + cols + " --out " + output_file,
            shell=True
        )
    else:
        p = call(
            plink + " --bfile " + input_file + " --score " + pgs_file + " " + cols + " --out " + output_file,
            shell=True
        )
    return 'Process is complete. Have a nice day!'


def download_pgs(
        *,
        prs_id,
):
    """
    Get PGS from pgscatalog.

    Parameters
    ----------
    prs_id : str
        the PRS ID in the pgscatalog.

    Returns
    -------
        a file containing the prs file

    """
    # make sure that the columns are present and matching
    resp = get('https://www.pgscatalog.org/rest/score/%s' % prs_id)
    prs_info = resp.json()
    if resp.status_code != 200 or not prs_info:
        raise Exception('The PRS ID might be wrong!')
    url = prs_info['ftp_scoring_file']
    prs_file = prs_id + '.gz'
    request.urlretrieve(url, prs_file)
    return prs_file
