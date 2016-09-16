import logging
import click
import os
import pdb
import sys

from .call import call as call_cmd
from .call import event_filter as event_filter_cmd
from .call import cluster_filter as cluster_filter_cmd


@click.group()
@click.option('--loglevel', default='INFO', help='level of logging')
@click.pass_context
def base(ctx, loglevel):
    setup_logging(loglevel)


def setup_logging(loglevel="INFO"):
    """
    Set up logging
    :param loglevel: loglevel to use, one of ERROR, WARNING, DEBUG, INFO (default INFO)
    :return:
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(levelname)s %(asctime)s %(funcName)s - %(message)s')
    logging.debug("Started log with loglevel %(loglevel)s" % {"loglevel": loglevel})


base.add_command(call_cmd)
base.add_command(event_filter_cmd)
base.add_command(cluster_filter_cmd)
