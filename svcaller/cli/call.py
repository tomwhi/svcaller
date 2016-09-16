import logging
import os

import click
import pysam

from svcaller.cli.calling.calling import event_filt, clust_filt
from svcaller.cli.calling.calling import DEL, INV, TRA, DUP


@click.command()
@click.option('--event-type', type=click.Choice([DEL, INV, TRA, DUP]), required=True)
@click.option('--output-name', type=str, default="SV_CallOutput", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def call(ctx, event_type, output_name, input_bam):
    logging.info("Running event calling: {}".format(event_type))

    logging.info("Applying event-specific filter...")
    event_filt_reads_bam_name = output_name + "_filt1.bam"
    event_filt_cmd(event_type, event_filt_reads_bam_name, input_bam)

    logging.info("Applying read-cluster filter...")
    clust_filt_reads_bam_name = output_name + "_filt2.bam"
    cluster_filt_cmd(clust_filt_reads_bam_name, event_filt_reads_bam_name)

    logging.info("Writing output...")
    with pysam.AlignmentFile(output_name + "_filtered_reads.bam", "wb", header=samfile.header) as outf:
        for read in cluster_filtered_reads:
            outf.write(read)

    #call_event(input_bam, output_name, event_type)


@click.command()
@click.option('--event-type', type=click.Choice([DEL, INV, TRA, DUP]), required=True)
@click.option('--output-bam', type=str, default="event_filtered_reads.bam", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def event_filter(ctx, event_type, output_bam, input_bam):
    event_filt_cmd(event_type, output_bam, input_bam)


def event_filt_cmd(event_type, output_bam, input_bam):
    logging.info("Running event filter: {}".format(event_type))

    samfile = pysam.AlignmentFile(input_bam, "rb")
    event_filtered_reads = event_filt(samfile, event_type)

    with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as outf:
        for read in event_filtered_reads:
            outf.write(read)


@click.command()
@click.option('--output-name', type=str, default="cluster_filtered_reads", required=False)
@click.argument('input-bam', type=click.Path(exists=True))
@click.pass_context
def cluster_filter(ctx, output_bam, input_bam):
    cluster_filter_cmd(output_bam, input_bam)


def cluster_filter_cmd(output_bam, input_bam):
    logging.info("Running cluster filter.")

    samfile = pysam.AlignmentFile(input_bam, "rb")
    reads = [r for r in list(samfile)]
    cluster_filtered_reads = clust_filt(reads, samfile, output_name)
    with pysam.AlignmentFile(output_bam, "wb", header=samfile.header) as outf:
        for read in cluster_filtered_reads:
            outf.write(read)