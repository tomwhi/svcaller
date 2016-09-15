import logging
import os

import click

from svcaller.cli.calling.calling import call_event
from svcaller.cli.calling.calling import DEL, INV, TRA, DUP


@click.command()
@click.option('--event-type', type=click.Choice([DEL, INV, TRA, DUP]), required=True)
@click.option('--output-name', type=str, default="SV_CallOutput", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def call(ctx, event_type, output_name, input_bam):
    logging.info("Running event calling: {}".format(event_type))

    # FIXME: Probably a more elegant way to do this:
    call_event(input_bam, output_name, event_type)
