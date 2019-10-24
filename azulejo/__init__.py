# -*- coding: utf-8 -*-
'''
azulejo -- tile phylogenetic space with subtrees
'''
import functools
import locale
import logging
import sys
from datetime import datetime
from pathlib import Path
#
# third-party imports
#
import click
import coverage
#
# package imports
#
from .version import version as VERSION
from .common import *
#
# start coverage
#
coverage.process_startup()
#
# set locale so grouping works
#
for localename in ['en_US', 'en_US.utf8', 'English_United_States']:
    try:
        locale.setlocale(locale.LC_ALL, localename)
        break
    except locale.Error:
        continue
#
# global constants
#
PROGRAM_NAME='azulejo'
AUTHOR = 'Joel Berendzen'
EMAIL = 'joelb@ncgr.org'
COPYRIGHT = """Copyright (C) 2019, NCGR. All rights reserved.
"""
LOG_DIR = 'logs'
LOG_PATH = Path('.') / LOG_DIR
# defaults for command line
DEFAULT_FILE_LOGLEVEL = logging.DEBUG
DEFAULT_STDERR_LOGLEVEL = logging.INFO
#
# Class definitions.
#
class CleanInfoFormatter(logging.Formatter):
    def __init__(self, fmt='%(levelname)s: %(message)s'):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        if record.levelno == logging.INFO:
            return record.getMessage()
        return logging.Formatter.format(self, record)
#
# time stamp for start
#
STARTTIME = datetime.now()
#
# global logger object
#
logger = logging.getLogger(PROGRAM_NAME)
#
# private context function
#
_ctx = click.get_current_context
#
# Helper functions
#
def init_dual_logger(file_log_level=DEFAULT_FILE_LOGLEVEL,
                     stderr_log_level=DEFAULT_STDERR_LOGLEVEL):
    '''Log to stderr and to a log file at different levels
    '''
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            global logger
            # find out the verbose/quiet level
            if _ctx().params['verbose']:
                _log_level = logging.DEBUG
            elif _ctx().params['quiet']:
                _log_level = logging.ERROR
            else:
                _log_level = stderr_log_level
            logger.setLevel(file_log_level)
            stderrHandler = logging.StreamHandler(sys.stderr)
            stderrFormatter = CleanInfoFormatter()
            stderrHandler.setFormatter(stderrFormatter)
            stderrHandler.setLevel(_log_level)
            logger.addHandler(stderrHandler)
            if _ctx().params['log']:  # start a log file in LOG_PATH
                logfile_path = LOG_PATH / (PROGRAM_NAME + '.log')
                if not LOG_PATH.is_dir():  # create LOG_PATH
                    try:
                        logfile_path.parent.mkdir(mode=0o755, parents=True)
                    except OSError:
                        logger.error('Unable to create log directory "%s"',
                                     logfile_path.parent)
                        raise OSError
                else:
                    if logfile_path.exists():
                        try:
                            logfile_path.unlink()
                        except OSError:
                            logger.error('Unable to remove log file "%s"',
                                         logfile_path)
                            raise OSError
                logfileHandler = logging.FileHandler(str(logfile_path))
                logfileFormatter = logging.Formatter(
                    '%(levelname)s: %(message)s')
                logfileHandler.setFormatter(logfileFormatter)
                logfileHandler.setLevel(file_log_level)
                logger.addHandler(logfileHandler)
            logger.debug('Command line: "%s"', ' '.join(sys.argv))
            logger.debug('%s version %s', PROGRAM_NAME, VERSION)
            logger.debug('Run started at %s', str(STARTTIME)[:-7])
            return f(*args, **kwargs)
        return wrapper
    return decorator


def init_user_context_obj(initial_obj=None):
    '''Put info from global options into user context dictionary
    '''
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            global config_obj
            if initial_obj is None:
                _ctx().obj = {}
            else:
                _ctx().obj = initial_obj
            ctx_dict = _ctx().obj
            if _ctx().params['verbose']:
                ctx_dict['logLevel'] = 'verbose'
            elif _ctx().params['quiet']:
                ctx_dict['logLevel'] = 'quiet'
            else:
                ctx_dict['logLevel'] = 'default'
            return f(*args, **kwargs)
        return wrapper
    return decorator


@click.group(epilog=AUTHOR + ' <' + EMAIL + '>. ' + COPYRIGHT)
@click.option('--warnings_as_errors', is_flag=True, show_default=True,
              default=False, help='Warnings cause exceptions.')
@click.option('-v', '--verbose', is_flag=True, show_default=True,
              default=False, help='Debug info to stderr.')
@click.option('-q', '--quiet', is_flag=True, show_default=True,
              default=False, help='Suppress logging to stderr.')
@click.option('--log/--no-log', is_flag=True, show_default=True,
              default=True, help='Write analysis in ./' + LOG_DIR + '.')
@click.version_option(version=VERSION, prog_name=PROGRAM_NAME)
@init_dual_logger()
@init_user_context_obj()
def cli(warnings_as_errors, verbose, quiet, log):
    """azulejo -- tiling gene in subtrees across phylogenetic space

    If COMMAND is present, and --no_log was not invoked,
    a log file named azulejo-COMMAND.log
    will be written in the ./logs/ directory.
    """
    if warnings_as_errors:
        logger.debug('Runtime warnings (e.g., from pandas) will cause exceptions')
        warnings.filterwarnings('error')


# other cli components import here
from .core import usearch_cluster, cluster_in_steps, clusters_to_histograms
from .analysis import analyze_clusters, compare_clusters

@cli.command()
def test():
    TESTDIR = '/home/localhost/joelb/preclust/'
    genomes = {'medtr': {'all': ['jemalong_A17', 'R108_HM340', 'HM004',
                                 'HM010', 'HM022', 'HM023',
                                 'HM034', 'HM050', 'HM056',
                                 'HM058', 'HM060', 'HM095',
                                 'HM125', 'HM129', 'HM185',
                                 'HM324',
                                 'jema+HM056',
                                 'jema+R108',
                                 'all'
                                ],
                         'refs': ['jemalong_A17']},
               'glyma':{'all': ['Williams',
                                'Lee',
                                'Zhonghuang',
                                'all',
                                #'Wm+glyso',
                                #'Wm+phavu',
                                #'Wm+jemalong',
                                ],
                        'ref': 'Williams'},
               'glyso': {'all': ['W05',
                                 'PI'],
                         'ref': 'W05'},
               'glyma+glyso': {'all': ['all'],
                               'ref': 'all'},
               'vigun': {'all': ['IT97K',
                                 'IT97K+CB5',
                                 'IT97K+Sanzi',
                                 'IT97K+Suvita',
                                 'IT97K+TZ30',
                                 'IT97K+UCR779',
                                 'IT97K+Xiabao',
                                 'IT97K+ZN016'
                                 ],
                         'ref': 'IT97K'}
               }
    species = 'vigun'
    #operation = 'compute'
    operation = 'single'
    operation = 'plot'
    #operation = 'condense'
    #operation = 'compare'
    print('doing species %s' %species)
    if operation == 'compute':
        for gnm in genomes[species]['all']:
            cluster_in_steps(TESTDIR + species + '/' + gnm,
                             'protein.faa'%(species, gnm),
                             16,
                             substrs='log/Substr.tsv',
                             dups='log/Dups.tsv',
                             min_id_freq=10)
    elif operation == 'single':
        gnm = genomes[species]['all'][0]
        usearch_cluster(TESTDIR + species + '/' +gnm +'/protein.faa',
                        0.984,
                        delete=False,
                        write_ids=True,
                        #do_calc=False,
                        subtrs='log/Substr.tsv',
                        dupss='log/Dups.tsv',
                        min_id_freq=10,
                        )
    elif operation == 'plot':
        namelist = []
        for gnm in genomes[species]['all']:
            namelist.append('%s/protein'%(gnm))
        analyze_clusters(TESTDIR + species,
                        namelist, species,
                        reference='%s/protein'%(
                                          genomes[species]['ref']),
                         on_id='vigun',
                         match_type='all')
    elif operation == 'condense':
        clusters_to_histograms(TESTDIR+'glyma+glyso/steven',
                               'steven_clusters.tsv')
    elif operation == 'compare':
        compare_clusters(TESTDIR+'glyma+glyso/steven/steven_clusters.tsv',
                         TESTDIR+'glyma+glyso/all/glyma+glyso.all-nr-984-ids.tsv')

if __name__ == '__main__':
    cli()